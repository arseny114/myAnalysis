import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt
import json
import os

# ==============================================================================
# 1. ЗАГРУЗКА ДАННЫХ
# ==============================================================================
print("Загрузка данных...")
df = pd.read_csv('ml_data.csv')
print(f"Загружено событий: {len(df)}")

# Разделяем на сигнал и фон
sig_df = df[df['is_signal'] == 1].copy()
bkg_df = df[df['is_signal'] == 0].copy()

print(f"Сигнал (сырых событий): {len(sig_df)}")
print(f"Фон (сырых событий): {len(bkg_df)}")

# ==============================================================================
# 2. НОРМАЛИЗАЦИЯ ВЕСОВ
# ==============================================================================
# Мы хотим, чтобы суммарный вес сигнала был равен суммарному весу фона.
# Это предотвратит доминирование огромного фона qq в функции потерь.

sum_w_sig = sig_df['weight'].sum()
sum_w_bkg = bkg_df['weight'].sum()

# Целевой суммарный вес (выбираем максимум из двух, или среднее, чтобы не терять статистику)
target_sum_w = max(sum_w_sig, sum_w_bkg) 

# Масштабируем веса
sig_df['scaled_weight'] = sig_df['weight'] * (target_sum_w / sum_w_sig)
bkg_df['scaled_weight'] = bkg_df['weight'] * (target_sum_w / sum_w_bkg)

# Опционально: умножим на 1000, чтобы веса были > 1 (XGBoost работает стабильнее с весами > 1)
scale_factor = 1000.0 / target_sum_w 
sig_df['scaled_weight'] *= scale_factor
bkg_df['scaled_weight'] *= scale_factor

print(f"Сумма весов сигнала после нормализации: {sig_df['scaled_weight'].sum():.2f}")
print(f"Сумма весов фона после нормализации: {bkg_df['scaled_weight'].sum():.2f}")

# Объединяем обратно
df_balanced = pd.concat([sig_df, bkg_df], ignore_index=True)

# ==============================================================================
# 3. ПОДГОТОВКА ПРИЗНАКОВ (FEATURES)
# ==============================================================================
# Исключаем мета-данные и целевую переменную
exclude_cols = ['process', 'is_signal', 'weight', 'scaled_weight']
features = [col for col in df.columns if col not in exclude_cols]

X = df_balanced[features].values
y = df_balanced['is_signal'].values
w = df_balanced['scaled_weight'].values

# Делим на Train и Test
X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
    X, y, w, test_size=0.3, random_state=42, stratify=y
)

print(f"Train set: {len(X_train)} событий")
print(f"Test set: {len(X_test)} событий")

# ==============================================================================
# 4. ОБУЧЕНИЕ XGBOOST
# ==============================================================================
print("\nНачало обучения XGBoost...")

# Гиперпараметры, оптимизированные для HEP (защита от переобучения)
params = {
    'objective': 'binary:logistic',
    'eval_metric': 'auc',
    'max_depth': 4,             # Небольшая глубина для физических переменных
    'learning_rate': 0.05,
    'subsample': 0.8,           # Сэмплирование строк
    'colsample_bytree': 0.8,    # Сэмплирование признаков
    'min_child_weight': 10,     # Учитывает scaled_weight
    'gamma': 0.1,
    'tree_method': 'hist',      # Быстрый алгоритм
    'device': 'cpu'
}

dtrain = xgb.DMatrix(X_train, label=y_train, weight=w_train, feature_names=features)
dtest = xgb.DMatrix(X_test, label=y_test, weight=w_test, feature_names=features)

model = xgb.train(
    params,
    dtrain,
    num_boost_round=500,
    evals=[(dtrain, 'train'), (dtest, 'eval')],
    early_stopping_rounds=30,
    verbose_eval=50
)

# Сохранение модели для дальнейшего использования в C++
model.save_model("xgb_zh_invisible.json")
print("\nМодель сохранена в xgb_zh_invisible.json")

# ==============================================================================
# 5. ОЦЕНКА КАЧЕСТВА И ГРАФИКИ
# ==============================================================================
y_pred_proba = model.predict(dtest)

auc_score = roc_auc_score(y_test, y_pred_proba, sample_weight=w_test)
print(f"ROC AUC на тестовой выборке: {auc_score:.4f}")

# Создаем директорию для сохранения графиков
os.makedirs("pdf_results", exist_ok=True)

# --- График ROC кривой ---
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba, sample_weight=w_test)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {auc_score:.3f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate (Background Efficiency)')
plt.ylabel('True Positive Rate (Signal Efficiency)')
plt.title('Receiver Operating Characteristic (Weighted)')
plt.legend(loc="lower right")
plt.grid(True, alpha=0.3)
plt.savefig("pdf_results/roc_curve.pdf")

# --- График важности признаков ---
xgb.plot_importance(model, max_num_features=15, importance_type='gain')
plt.title("Feature Importance (Gain)")
plt.tight_layout()
plt.savefig("pdf_results/feature_importance.pdf")

# --- Распределение скоринга (BDT Output) ---
plt.figure(figsize=(8, 6))
# Сигнал
plt.hist(y_pred_proba[y_test == 1], bins=50, weights=w_test[y_test == 1], 
         alpha=0.6, color='red', label='Signal (qqHinvi)', density=True)
# Фон
plt.hist(y_pred_proba[y_test == 0], bins=50, weights=w_test[y_test == 0], 
         alpha=0.6, color='blue', label='Background', density=True)

plt.xlabel('BDT Score')
plt.ylabel('Weighted Density')
plt.title('BDT Output Distribution')
plt.legend()
plt.yscale('log')
plt.grid(True, alpha=0.3)
plt.savefig("pdf_results/bdt_output_dist.pdf")
