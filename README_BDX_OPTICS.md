# Ottica dei cristalli BGO BDX

Repository: [marcospreafico/clas12Tags](https://github.com/marcospreafico/clas12Tags/tree/optical_photons_simulation),
ramo **`optical_photons_simulation`**, derivato dal ramo `frederic` al commit
`1d2313f8cc701e43bf461bee44ea12f43c7d402d`. Per clonarlo:

```bash
git clone --branch optical_photons_simulation --single-branch \
  https://github.com/marcospreafico/clas12Tags.git clas12Tags
```

Simulazione GEMC dei singoli cristalli BGO, con otto variation `bgo1` ... `bgo8`,
estratte da `experiments/bdx/geometry_crs60.pl` del ramo `frederic`.
I vertici `G4GenericTrap` originali sono conservati. Tutti i tipi sono lunghi
**240 mm**, lungo z: base minore a -120 mm, base maggiore a +120 mm.
Nessuna geometria o proprietà dei materiali di MINE è importata.

La produzione registra il **primo arrivo di ogni fotone ottico su tutta la base
maggiore**, con coordinate, direzione, tempo e storia. Il collettore `bdx_flux`
è una sottile estrusione della stessa base trapezoidale, con ingresso adattato
all'indice BGO: il fotone viene contato e fermato, senza applicare PDE o tagli 3×3.
Il resto del cristallo è rivestito con il modello di vernice bianca.

## Produzione

Richiede Docker Desktop attivo e Python >= 3.9. Matplotlib serve solo all'analisi
(`python3 -m pip install -r requirements-bdx-optics.txt`). L'immagine GEMC è fissata
per digest. `--native` usa invece un ambiente GEMC/HIPO locale già configurato.

```bash
cd /Users/simone/Work/JLAB/BDX/bdx-crystals/clas12Tags
bash scripts/build_bdx_optics.sh

python3 scripts/run_bdx_optics.py \
  --types 1 2 3 4 5 6 7 8 --distances-mm 20 80 140 200 --events 10 \
  --set light_yield_per_MeV=800
```

Non viene lanciata alcuna analisi o prodotto alcun grafico. Ogni esecuzione crea
una nuova directory in `build/bdx_optics/`; `--output` permette di sceglierla,
sempre all'interno del clone. Directory esistenti sono rifiutate.
`--dry-run` prepara gli input senza simulare. Le otto geometrie vengono rigenerate
negli input della campagna, senza modificare il Perl originale.

Per ogni punto `bgoN/d...mm/` vengono scritti:

- `optical.hipo`: dati grezzi, con le banche descritte in [OUTPUT_BANKS.md](experiments/bdx_optics/OUTPUT_BANKS.md);
- `events.csv`: riepilogo per evento, duplicato leggero per controllo;
- `banks.json`: inventario delle banche e numero di righe, dopo rilettura e verifica;
- `gemc.log`, `hipo_check.log`, `run.gcard`, `optics.conf`: log e input del punto.

`campaign.json` registra stato, semi distinti, configurazione, geometria, immagine,
hash di codice, eseguibile, CSV e HIPO. `inputs/` conserva il riferimento originale.
Un errore di trasporto, bilancio o integrità HIPO interrompe la campagna.
Sono salvati anche gli eventi con zero arrivi. Non esistono tagli sul tempo,
numero di passi o numero di riflessioni dei fotoni.

Default: muoni negativi da 4 GeV/c lungo +x, partenza (-100,0,z) mm e
`z=120 mm-distanza`. `--particle` e `--momentum-GeV` sono configurabili.
La distanza della scansione è assiale, dalla base maggiore; il cammino del
fotone viene salvato separatamente. Tutti i secondari della particella contribuiscono
alla deposizione e alla generazione ottica nel BGO.

## Analisi indipendente

Sostituire `NOME_CAMPAGNA` con la directory stampata al termine della produzione.

```bash
# Raccolta su tutta la base maggiore:
python3 scripts/analyze_bdx_optics.py build/bdx_optics/NOME_CAMPAGNA

# Sugli stessi fotoni: quadrato centrale 3×3 mm²:
python3 scripts/analyze_bdx_optics.py build/bdx_optics/NOME_CAMPAGNA --square-mm 3
```

Ogni analisi crea una directory nuova **fuori dalla produzione**:
`build/bdx_optics/NOME_CAMPAGNA_analysis/SELEZIONE_TIMESTAMP/`.
`--output PERCORSO_NUOVO` la sceglie esplicitamente. I dati grezzi non sono modificati.
Sono prodotti `collection.csv`, `analysis.json`, `collection.svg`, `collection.png`
e una copia CSV delle banche lette in `banks/`. Matplotlib è richiesto qui.

Il quadrato è centrato nel baricentro di area del trapezio, con assi x/y;
la selezione usa `abs(lx_mm)<=1.5` e `abs(ly_mm)<=1.5`, solo `origin=1`
(scintillazione). Nessuna PDE è applicata. Per definire analisi differenti,
esportare le banche senza eseguire questo esempio di analisi:

```bash
python3 scripts/export_bdx_optics.py \
  build/bdx_optics/NOME_CAMPAGNA/bgo1/d20mm/optical.hipo \
  --output build/bdx_optics/NOME_EXPORT
```

L'export contiene `Config.csv`, `Event.csv`, `Edep.csv`, `Flux.csv` e `banks.json`.
`record_index` unisce le tabelle evento per evento. Le banche HIPO restano
leggibili direttamente da normali strumenti HIPO; lo schema è nel file.

Normalizzazioni dell'analisi d'esempio, con somme su **tutti** gli eventi:

- raccolta geometrica: `somma(arrivi selezionati scint) / somma(fotoni generati scint)`;
- fotoni per MeV: `somma(arrivi selezionati scint) / somma(Edep_MeV)`;
- confronto con la resa impostata: il precedente rapporto diviso per `light_yield_per_MeV`.

Si salva anche `somma(generati)/somma(Edep)`. Il primo rapporto usa i fotoni
realmente generati, il terzo la resa impostata: differiscono per le fluttuazioni
statistiche. L'incertezza principale è il SEM del rapporto sugli eventi indipendenti;
un solo evento non fornisce questo errore. L'intervallo Wilson 68% descrive soltanto
la precisione del conteggio dei fotoni. Gli errori del modello non sono inclusi.
Abbassare la resa riduce il tempo di calcolo e aumenta le fluttuazioni, senza pesi nascosti.

**Interpretazione della selezione a posteriori:** la base intera è un assorbitore
ideale al primo arrivo. Un taglio 3×3 descrive la frazione dei primi arrivi nel
quadrato. Un vero piccolo SiPM con il resto della base verniciato può raccogliere
anche fotoni prima riflessi fuori dal quadrato: è una diversa condizione al contorno.
Le vecchie produzioni con SiPM 3×3 e assorbimento nullo non contengono i primi
arrivi sulla base intera; occorre una nuova produzione. L'analisi legge ancora i
vecchi `events.csv` ma rifiuta di applicarvi nuovi tagli spaziali.

## Parametri ottici

File `experiments/bdx_optics/geometry_only.conf`; gli override `--set` vengono
salvati negli input e in `BDX::Config`. Le nuove produzioni usano **10 m** di
lunghezza di assorbimento BGO, valore comunicato dall'utente il 24 novembre 2026.

| Parametro | Default | Significato |
|---|---:|---|
| `absorption_length_mm` | **10000** | Lunghezza di assorbimento nel volume BGO; 0 lo disattiva esplicitamente |
| `wrap_reflectivity` | 0.98 | Probabilità assoluta di riflessione a ogni incontro con la vernice |
| `diffuse_length_mm` | 40 | Zona Lambertiana opaca dalla base maggiore |
| `transition_length_mm` | 100 | Miscela lineare da Lambertiana a speculare |
| `specular_length_mm` | 100 | Zona speculare fino alla base minore |
| `light_yield_per_MeV` | 8500 | Resa proporzionale di riferimento |
| `refractive_index` | 2.15 | Indice costante, collettore adattato |
| `decay_time_ns` | 300 | Decadimento esponenziale |
| `wavelength_nm` | 480 | Picco stretto equivalente; coefficienti ottici costanti |
| `cherenkov` | 0 | Con 1 produzione e arrivi Cherenkov restano separati |

Per `d=120 mm-z` sulla superficie verniciata, `f_spec=0` fino a 40 mm,
`(d-40)/100` fra 40 e 140 mm, `1` oltre 140 mm. Le probabilità sono
`P_ass=1-R`, `P_spec=R*f_spec`, `P_diff=R*(1-f_spec)`: R non cambia con la miscela.
La vernice è opaca. Le tre lunghezze devono sommare alla lunghezza fisica 240 mm.
Le facce laterali e la base minore sono verniciate; tutta la base maggiore è aperta
verso il collettore. L'assorbimento nel BGO segue `exp(-cammino/10000 mm)`.

Si trascurano accoppiamento ottico, package/PDE del SiPM, rumore, saturazione,
dispersione cromatica, Rayleigh/Mie/WLS e non proporzionalità della scintillazione.
Il picco stretto non è uno spettro BGO misurato. Nessun materiale ottico MINE è usato.
Riferimenti e provenienza: [reference/bdx_optics/README.md](reference/bdx_optics/README.md).

## Verifiche numeriche

```bash
bash scripts/test_bdx_optics.sh
```

I test confrontano gli otto solidi con il Perl originale, copertura esatta del
collettore, trasporto analitico, riflettività e miscela angolare, TIR, assorbimento
esponenziale, propagazione oltre 100 passi e integrità HIPO anche oltre 8 MiB.
I test ideali senza assorbimento lo disattivano esplicitamente.

Geant4 10.7.4 richiede particolare attenzione alle facce non piane di
`G4GenericTrap`. Il solido conserva vertici, superfici e normali originali, usando
limiti conservativi per le distanze senza direzione e radici quadratiche stabili
per le intersezioni. Dopo una riflessione, la posizione viene spostata strettamente
all'interno di circa 1 nm e sincronizzata con il volume corrente. La correzione
massima ammessa è 1 micrometro più `boundary_push_mm`; oltre tale limite il run
fallisce. Lo stesso controllo ripara arrotondamenti dei fotoni nati sul bordo,
senza ricampionare l'emissione nativa. Contatori e massimi sono salvati in
`BDX::Event` e nel CSV. La correzione non cambia energia o probabilità fisiche.

Il modello usa il `G4RunManager` sequenziale di questo GEMC. Scansioni parallele
richiedono processi, directory e semi indipendenti. Versione Geant4 verificata: 10.7.4.
Le gcards vanno usate con il lanciatore: la configurazione `BDX_OPTICS_CONFIG`
attiva il codice ottico dedicato, che non è attivo per altri detector GEMC.

## Farm JLab e molti file HIPO

I comandi di compilazione nativa, pianificazione e sottomissione SWIF2 sono in
[farm/README.md](farm/README.md). Ogni job contiene un solo tipo e una sola
distanza, con numero di eventi configurabile. `scripts/analyze_bdx_hipos.py`
aggrega campagne, glob o liste di HIPO, con selezione opzionale `--square-mm 3`.
