import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import pyreadr
import os

st.set_page_config(page_title="Symulacja: Regresja Porządkowa", layout="wide")

st.title("Symulacja Monte Carlo – Regresja Porządkowa")

st.sidebar.header("Parametry Symulacji")
n_reps = st.sidebar.number_input("Liczba powtórzeń (MC):", min_value=5, max_value=200, value=20, step=5)
n = st.sidebar.number_input("Próba ucząca (n):", min_value=50, max_value=1000, value=100, step=50)
m = st.sidebar.number_input("Próba testowa (m):", min_value=100, max_value=2000, value=500, step=100)
p = st.sidebar.number_input("Liczba zmiennych (p):", min_value=5, max_value=200, value=30, step=5)

model_type = st.sidebar.selectbox(
    "Mechanizm generowania (Y*):",
    options=["M1", "M2", "M3"],
    format_func=lambda x: {
        "M1": "Liniowy z bł. normalnym (M1)",
        "M2": "Liniowy z bł. Cauchy'ego (M2)",
        "M3": "Nieliniowy wykładniczy (M3)"
    }[x]
)

beta_type = st.sidebar.selectbox(
    "Wektor współczynników:",
    options=["betaA", "betaB"],
    format_func=lambda x: "beta_A (3 zmienne)" if x == "betaA" else "beta_B (10 zmiennych)"
)

rho = st.sidebar.slider("Korelacja (rho):", min_value=0.0, max_value=0.9, value=0.5, step=0.1)

run_button = st.sidebar.button("Uruchom symulację", type="primary")

if run_button:
    output_file = "temp_results.RData"
    
    cmd = [
        "Rscript", "sim.R",
        str(n), str(m), str(p), str(rho), str(n_reps),
        model_type, beta_type, output_file
    ]
    
    with st.spinner("Uruchamianie silnika R i wykonywanie symulacji Monte Carlo..."):
        try:
            # Przechwytywanie komunikatów wyjścia oraz błędów
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            results = pyreadr.read_r(output_file)
            df_pred = results["boxplot_prediction"]
            df_est = results["boxplot_estimation"]
            
            if os.path.exists(output_file):
                os.remove(output_file)

            st.success("Symulacja zakończona!")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Błąd Predykcji")
                fig1, ax1 = plt.subplots(figsize=(6, 4))
                df_pred.boxplot(column='error', by='model', ax=ax1, grid=False)
                plt.title("")
                plt.suptitle("")
                plt.xlabel("Model")
                plt.ylabel("Błąd (1 - ranking acc)")
                st.pyplot(fig1)

            with col2:
                st.subheader("Błąd Estymacji")
                fig2, ax2 = plt.subplots(figsize=(6, 4))
                df_est.boxplot(column='error', by='model', ax=ax2, grid=False)
                plt.title("")
                plt.suptitle("")
                plt.xlabel("Model")
                plt.ylabel("Błąd sferyczny ||b_n - b_true||")
                st.pyplot(fig2)

            st.markdown("---")
            st.subheader("Podsumowanie Statystyczne")
            
            sum_pred = df_pred.groupby('model')['error'].agg(['mean', 'std']).reset_index()
            sum_pred.columns = ['Model', 'Predykcja (Średnia)', 'Predykcja (Odch. Std.)']
            
            sum_est = df_est.groupby('model')['error'].agg(['mean', 'std']).reset_index()
            sum_est.columns = ['Model', 'Estymacja (Średnia)', 'Estymacja (Odch. Std.)']
            
            summary_table = pd.merge(sum_pred, sum_est, on='Model')
            st.dataframe(summary_table, use_container_width=True)

        except subprocess.CalledProcessError as e:
            st.error("Wystąpił błąd podczas wykonywania skryptu R:")
            st.code(e.stderr if e.stderr else e.stdout, language="R")
else:
    st.info("Ustaw parametry w panelu bocznym i kliknij **Uruchom symulację**.")
