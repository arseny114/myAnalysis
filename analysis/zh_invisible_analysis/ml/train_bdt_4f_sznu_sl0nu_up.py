import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt
import json
import os

# Настройка стиля графиков
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 12

print("=" * 80)
print("1. ЗАГРУЗКА И ФИЛЬТРАЦИЯ ДАННЫХ")
print("=" * 80)
print("Загрузка данных из ml_data.csv...")
df = pd.read_csv('ml_data.csv')
print(f"Всего событий в CSV: {len(df)}")

# Фильтруем только сигнал и целевой фон
signal_process = 'E240_qqHinvi'
background_process = 'E240_4f_sznu_sl0nu_up'

sig_df = df[df['process'] == signal_process].copy()
bkg_df = df[df['process'] == background_process].copy()

print(f"Сигнал ({signal_process}): {len(sig_df)} событий")
print(f"Фон ({background_process}): {len(bkg_df)} событий")

if len(sig_df) == 0 or len(bkg_df) == 0:
    raise ValueError("Один из процессов отсутствует в ml_data.csv. Убедитесь, что вы запустили main.cpp с флагом --export-csv и указали оба ROOT-файла.")

print("\n" + "=" * 80)
print("2. БАЛАНСИРОВКА 50/50 (БЕЗ ВЕСОВ)")
print("=" * 80)
# Берем минимальное количество событий из двух классов для строгого баланса 50/50
min_events = min(len(sig_df), len(bkg_df))
print(f"Балансировка: ограничиваем каждый класс до {min_events} событий (random_state=42)")

sig_df_balanced = sig_df.sample(n=min_events, random_state=42)
bkg_df_balanced = bkg_df.sample(n=min_events, random_state=42)

df_balanced = pd.concat([sig_df_balanced, bkg_df_balanced], ignore_index=True)
print(f"Итоговый размер сбалансированного датасета: {len(df_balanced)} событий")

print("\n" + "=" * 80)
print("3. ПОДГОТОВКА ПРИЗНАКОВ (FEATURES)")
print("=" * 80)
# Порядок должен строго совпадать с заполнением в main.cpp и списком в zh_invisible_analysis.h
features = [
    "invMass", "cosThetaZ", "deltaR", "cosTheta1", "cosTheta2", 
    "jet1_E", "jet2_E", "met_jet", "dijetEnergy", "deltaTheta", 
    "deltaPhi", "pmiss_mag", "cosThetaPmiss"
]

X = df_balanced[features].values
y = df_balanced['is_signal'].values

# Веса не используем при обучении (равны 1), но сохраним оригинальные веса для оценки физического качества на тесте
w_physical = df_balanced['weight'].values 

X_train, X_test, y_train, y_test, w_train_physical, w_test_physical = train_test_split(
    X, y, w_physical, test_size=0.3, random_state=42, stratify=y
)

print(f"Train set: {len(X_train)} событий")
print(f"Test set: {len(X_test)} событий")

print("\n" + "=" * 80)
print("4. ОБУЧЕНИЕ XGBOOST (БЕЗ ВЕСОВ)")
print("=" * 80)
params = {
    'objective': 'binary:logistic',
    'eval_metric': 'auc',
    'max_depth': 5,
    'learning_rate': 0.05,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'min_child_weight': 5,
    'gamma': 0.1,
    'tree_method': 'hist',
    'device': 'cpu'
}

dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=features)
dtest = xgb.DMatrix(X_test, label=y_test, feature_names=features)

model = xgb.train(
    params,
    dtrain,
    num_boost_round=1000,
    evals=[(dtrain, 'train'), (dtest, 'eval')],
    early_stopping_rounds=50,
    verbose_eval=100
)

model_path = "xgb_zh_4f_sznu_sl0nu_up.json"
model.save_model(model_path)
print(f"\nМодель сохранена в {model_path}")

print("\n" + "=" * 80)
print("5. ОЦЕНКА КАЧЕСТВА И ГРАФИКИ")
print("=" * 80)
y_pred_proba = model.predict(dtest)

# Оцениваем AUC двумя способами:
# 1. Математический (без весов, просто качество разделения выборки 50/50)
auc_unweighted = roc_auc_score(y_test, y_pred_proba)
# 2. Физический (с весами, чтобы понять, как это работает на реальных сечениях)
auc_weighted = roc_auc_score(y_test, y_pred_proba, sample_weight=w_test_physical)

print(f"ROC AUC (Unweighted, 50/50 test): {auc_unweighted:.4f}")
print(f"ROC AUC (Weighted, physical test): {auc_weighted:.4f}")

os.makedirs("pdf_results_4f_sznu_sl0nu_up", exist_ok=True)

# --- График ROC кривой ---
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {auc_unweighted:.3f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate (4f_sznu_sl0nu_up Efficiency)')
plt.ylabel('True Positive Rate (qqHinvi Efficiency)')
plt.title('ROC Curve: qqHinvi vs 4f_sznu_sl0nu_up (Unweighted Training)')
plt.legend(loc="lower right")
plt.savefig("pdf_results_4f_sznu_sl0nu_up/roc_curve.pdf")
plt.close()

# --- График важности признаков ---
xgb.plot_importance(model, max_num_features=14, importance_type='gain', height=0.8)
plt.title("Feature Importance (Gain) - qqHinvi vs 4f_sznu_sl0nu_up")
plt.tight_layout()
plt.savefig("pdf_results_4f_sznu_sl0nu_up/feature_importance.pdf")
plt.close()

# --- Распределение скоринга (BDT Output) ---
plt.figure(figsize=(8, 6))
# Сигнал
plt.hist(y_pred_proba[y_test == 1], bins=50, 
         alpha=0.6, color='red', label='Signal (qqHinvi)', density=True)
# Фон
plt.hist(y_pred_proba[y_test == 0], bins=50, 
         alpha=0.6, color='blue', label='Background (4f_sznu_sl0nu_up)', density=True)
plt.xlabel('BDT Score')
plt.ylabel('Density')
plt.title('BDT Output Distribution (Test Set)')
plt.legend()
plt.yscale('log')
plt.savefig("pdf_results_4f_sznu_sl0nu_up/bdt_output_dist.pdf")
plt.close()

print("Графики сохранены в директорию pdf_results_4f_sznu_sl0nu_up/")
