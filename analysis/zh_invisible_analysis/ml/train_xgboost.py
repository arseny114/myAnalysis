#!/usr/bin/env python3
# train_xgboost.py
import argparse
import os
import numpy as np
import pandas as pd
import uproot
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, classification_report
import m2cgen as m2c

# Порядок фич должен строго совпадать с zh_invisible_analysis.h::MLFeatures
FEATURES = [
    "invMass", "recoilMass", "cosThetaZ", "deltaR",
    "ej1", "ej2", "eta_j1", "eta_j2",
    "pt_jj", "met_pfo", "pmag_miss", "costh_miss",
    "energy_asym", "nconst_j1", "nconst_j2"
]

def main():
    parser = argparse.ArgumentParser(description="Обучение XGBoost для ZH -> invisible")
    parser.add_argument("-s", "--signal", required=True, help="ROOT файл с сигналом (H->inv)")
    parser.add_argument("-b", "--background", required=True, help="ROOT файл с фоном (H->inclusive/other)")
    parser.add_argument("-o", "--output", default="../ml/xgboost_model.h", help="Путь к выходному C++ хедеру")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print("📦 Загрузка данных...")
    df_s = uproot.open(args.signal)["ml_features"].arrays(FEATURES + ["label"], library="pd")
    df_s["label"] = 1
    df_b = uproot.open(args.background)["ml_features"].arrays(FEATURES + ["label"], library="pd")
    df_b["label"] = 0
    df = pd.concat([df_s, df_b], ignore_index=True)
    print(f"  Всего событий: {len(df)} (Сигнал: {len(df[df.label==1])}, Фон: {len(df[df.label==0])})")

    X = df[FEATURES].values
    y = df["label"].values
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    print("🧠 Обучение XGBoost...")
    model = xgb.XGBClassifier(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, eval_metric="logloss",
        random_state=42, n_jobs=-1
    )
    model.fit(X_train, y_train, eval_set=[(X_test, y_test)], verbose=False)

    # Оценка
    y_proba = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_proba)
    print(f"\n📊 ROC AUC на тестовой выборке: {auc:.4f}")
    print(classification_report(y_test, y_proba > 0.5))

    # Важность фич
    print("\n🔍 Feature Importances:")
    for feat, imp in zip(FEATURES, model.feature_importances_):
        print(f"  {feat:15s}: {imp:.4f}")

    # Экспорт в C++ хедер
    print("\n💾 Генерация C++ хедера...")
    code = m2c.export_to_c(model)
    
    # m2cgen генерирует void score(...), оборачиваем под наш интерфейс
    header = """#ifndef XGBOOST_MODEL_H
#define XGBOOST_MODEL_H

#include <cmath>
#include <vector>

""" + code.replace("void score", "static void xgb_model_score") + """

inline double predict_xgboost(const double* features) {
    double proba[2];
    xgb_model_score(features, proba);
    // proba[0] = P(фон), proba[1] = P(сигнал)
    return proba[1];
}

#endif // XGBOOST_MODEL_H
"""
    with open(args.output, "w") as f:
        f.write(header)
    print(f"✅ Модель сохранена в: {args.output}")

if __name__ == "__main__":
    main()
