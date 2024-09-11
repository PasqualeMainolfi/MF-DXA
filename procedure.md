# MF-DXA Multifractal Detrended Cross-Correlation Analysis

1. Preparazione dei dati:
   - Inizia con due serie temporali x(i) e y(i), dove i = 1, 2, ..., M.
   - Assicurati che entrambe le serie abbiano la stessa lunghezza M.

2. Calcolo dei profili:
   - Calcola i profili X(i) e Y(i) sottraendo le medie dalle serie originali e sommando cumulativamente:
     X(i) = Σ[k=1 to i] [x(k) - x']
     Y(i) = Σ[k=1 to i] [y(k) - y']
   - Dove x' e y' sono le medie di x e y rispettivamente.

3. Divisione in segmenti:
   - Dividi entrambe le serie in Ns = int(N/s) segmenti non sovrapposti di lunghezza s.

4. Calcolo della covarianza locale:
   - Per ogni segmento v, calcola la covarianza locale F²xy(s,v) utilizzando polinomi di ordine m (solitamente lineari o quadratici) per detrend entrambe le serie.
   - La formula generale è:
     F²xy(s,v) = (1/s) Σ[k=1 to s] {[X((v-1)s+k) - X̃v(k)] * [Y((v-1)s+k) - Ỹv(k)]}
   - Dove X̃v(k) e Ỹv(k) sono i trend locali stimati.

5. Calcolo della funzione di fluttuazione:
   - Calcola la funzione di fluttuazione q-esima:
     Fq(s) = {1/Ns Σ[v=1 to Ns] sign[F²xy(s,v)] |F²xy(s,v)|^(q/2)}^(1/q)
   - Ripeti questo calcolo per diversi valori di q (sia positivi che negativi).

6. Ripetizione per diverse scale:
   - Ripeti i passaggi 3-5 per diverse lunghezze di segmento s.

7. Analisi della legge di scala:
   - Rappresenta graficamente log(Fq(s)) vs log(s) per ogni valore di q.
   - Se le serie sono cross-correlate in modo multiscala, dovresti osservare una relazione lineare.
   - Il coefficiente angolare di questa relazione lineare è l'esponente di Hurst generalizzato h(q).

8. Calcolo dello spettro multiscala:
   - Calcola lo spettro multiscala τ(q) = qh(q) - 1
   - Calcola la dimensione frattale generalizzata D(q) = τ(q) / (q-1)

9. Interpretazione dei risultati:
   - Analizza come h(q), τ(q) e D(q) variano con q per caratterizzare le proprietà multiscala della cross-correlazione tra le due serie.

Note:

- Per serie temporali molto lunghe, potrebbe essere necessario utilizzare tecniche di ottimizzazione computazionale.
- La scelta dell'ordine del polinomio per il detrending può influenzare i risultati, quindi potrebbe essere necessario sperimentare con diversi ordini.
- È importante verificare la robustezza dei risultati variando i parametri dell'analisi, come il range di s e q.