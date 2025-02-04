# 📌 Analisi delle Performance di un Modello Data-Driven per la Suscettibilità ai Crolli in Roccia

## 📖 Descrizione
Questo repository contiene il codice per l'analisi delle performance di un modello data-driven per la suscettibilità ai crolli in roccia nella Regione Valle d'Aosta. Il progetto confronta le performance del modello Generalized Additive Model (GAM) con la metodologia Multivariate Environmental Similarity Surface (MESS), valutando l'effetto della dissimilarità tra le aree di training e test.

L'analisi è suddivisa in due principali aree di studio:

Area 1: Comunità Montane Mont Emilius e Mont Cervin.
Area 2: Comunità Montane Evançon, Monte Rosa e Walser Alta Valle del Lys.
📄 Il lavoro si basa sulla mia tesi magistrale in Scienze della Terra, con focus sull’integrazione di variabili geologiche, topografiche e meteoclimatiche per migliorare la previsione del rischio di crolli in roccia.

## 📂 Struttura del Repository
bash
Copia
Modifica
├── data/                     # Dati di input utilizzati per il modello
│   ├── Area1_Rmodel.csv
│   ├── Area2_Rmodel.csv
│
├── scripts/                  # Script R per l'analisi
│   ├── MESS_comparison.R     # Script principale per l'analisi MESS e GAM
│
├── plots/                    # Output grafici generati dall'analisi
│
├── results/                  # Risultati finali e file CSV generati
│   ├── mess_auroc_CV.csv     # Risultati del cross-validation MESS-AUROC
│
├── Tesi_VdA.pdf              # Documento completo della tesi
│
└── README.md                 # Questo file
## 🛠 Installazione
1️⃣ Prerequisiti
R (≥ 4.0)
Pacchetti R: ggplot2, raster, sf, tidyverse, mgcv, dismo, sperrorest, pROC
2️⃣ Installazione dei Pacchetti
Apri RStudio e installa i pacchetti richiesti:

r
Copia
Modifica
install.packages(c("ggplot2", "raster", "sf", "tidyverse", "mgcv", "dismo", "sperrorest", "pROC"))
🚀 Esecuzione dell'Analisi
1️⃣ Carica lo script principale in R:

r
Copia
Modifica
source("scripts/MESS_comparison.R")
2️⃣ Assicurati che i file CSV siano presenti nella cartella data/.

3️⃣ Lo script esegue le seguenti operazioni:

Caricamento dei dati (Area1_Rmodel.csv e Area2_Rmodel.csv)
Trasformazione delle variabili con cutoff specifici
Modellazione della suscettibilità ai crolli con GAM
Calcolo delle metriche MESS e AUROC
Validazione con cross-validation random e K-means
Generazione di grafici per l’interpretazione dei risultati
## 📊 Output e Risultati
I risultati sono salvati nella cartella results/ e includono:

mess_auroc_CV.csv → Performance del modello in cross-validation.
Grafici nella cartella plots/, inclusi:
CV_mess_vs_diff_perf.jpg
boxplot_auroc_scores.jpg
A1_kmean_divided_train.jpg
## 📌 Possibili Estensioni
Espansione del dataset: aggiungere nuove aree geografiche o aggiornare i dati meteoclimatici.
Testare altri algoritmi di Machine Learning: come Random Forest o XGBoost per confrontare le prestazioni con GAM.
Automatizzare l'analisi con un workflow completo in RMarkdown o Jupyter Notebook.
## 📜 Licenza
Questo progetto è distribuito sotto licenza MIT - puoi liberamente modificarlo e condividerlo.
