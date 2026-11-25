# Banche della produzione ottica BDX — schema 2

File per punto: `bgoN/d...mm/optical.hipo`. Gli schemi sono incorporati nel file
HIPO e definiti in `source/src/BdxOpticalSchema.h`. Le quattro banche dedicate
usano gruppo 4601, item 1–4. `I` = intero 32 bit, `D` = double 64 bit.

| Banca | Righe per evento | Contenuto |
|---|---:|---|
| `BDX::Config` | 1 | Parametri ottici, geometria, tipo, seme, versione dello schema |
| `BDX::Event` | 1, anche senza arrivi | Energia depositata, fotoni generati, raccolti e persi |
| `BDX::Edep` | Una per passo sorgente con Edep>0 o fotoni prodotti | Deposizioni della particella/secondari e conteggi ottici associati |
| `BDX::Flux` | Una per primo arrivo ottico sulla base maggiore | Coordinate, direzione, energia/tempo e storia del singolo fotone |

Le righe nelle banche dello stesso evento HIPO si riferiscono allo stesso evento
fisico. `event_index` è l'ID Geant4 da zero; `event` corrisponde a
`RUN::config.event`, normalmente da uno. L'esportatore CSV aggiunge `record_index`
a **tutte** le tabelle, progressivo da zero nel singolo file, per unirle senza
ambiguità. Quando si uniscono più file, aggiungere anche l'identità del file/punto.
`source_step` unisce una riga Flux alla riga Edep dello stesso evento.

## `BDX::Config` — item 1

Tutti i parametri sono ripetuti per evento, senza dipendere da file esterni.

| Campi | Tipo | Significato/unità |
|---|---|---|
| `schema_version` | I | 2 |
| `crystal_type` | I | 1–8 nelle produzioni; 0 riservato a test |
| `seed` | I | Seme iniziale del punto di scansione, non un nuovo seme a ogni evento |
| `collector` | I | 1 = intera base maggiore al primo arrivo; 0 = apertura dei test analitici |
| `wrap_model` | I | 1 = vernice con miscela spaziale; 2 = speculare di confronto; 0 = superficie senza vernice |
| `cherenkov` | I | 0/1, processo disattivato/attivato |
| `light_yield_per_MeV` | D | Resa impostata, fotoni/MeV |
| `refractive_index` | D | Indice costante BGO e ingresso collettore |
| `absorption_length_mm` | D | **10000 mm** di default; 0 disattiva esplicitamente l'assorbimento nel BGO |
| `decay_time_ns`, `wavelength_nm` | D | Decadimento della scintillazione e picco stretto equivalente |
| `wrap_reflectivity` | D | R assoluta a ogni incontro, default 0.98 |
| `diffuse_length_mm`, `transition_length_mm`, `specular_length_mm` | D | Lunghezze delle tre zone, default 40/100/100 mm dalla base maggiore |
| `boundary_push_mm` | D | Correzione numerica interna di riferimento, default 1e-6 mm |
| `half_length_mm` | D | Semilunghezza lungo z, 120 mm |
| `Bg_mm`, `Bm_mm`, `H_mm` | D | Semibasi a y=+H/-H e semialtezza della faccia maggiore |
| `bg_mm`, `bm_mm`, `h_mm` | D | Corrispondenti semidimensioni della faccia minore |
| `center_y_mm` | D | Baricentro della base maggiore, `H*(Bg-Bm)/(3*(Bg+Bm))` |
| `plane_z_mm` | D | Piano di raccolta, +120 mm |

I vertici sono ricostruibili dalle semidimensioni: a -L sono
`(-bm,-h),(-bg,h),(bg,h),(bm,-h)`; a +L
`(-Bm,-H),(-Bg,H),(Bg,H),(Bm,-H)`. Il manifest della campagna conserva anche i
vertici espliciti e gli hash della sorgente originale e del codice.

## `BDX::Event` — item 2

| Campi | Tipo | Significato/unità |
|---|---|---|
| `event`, `event_index` | I | Numeri evento GEMC e Geant4 |
| `primary_z_mm`, `nominal_distance_mm` | D | z iniziale del fascio e distanza assiale `120-z` |
| `edep_MeV` | D | Edep nel BGO da tutte le tracce **non ottiche**, inclusi i secondari |
| `edep_axial_distance_mm` | D | Distanza assiale media, pesata con Edep, usando il punto medio di ciascun passo |
| `edep_euclidean_distance_mm` | D | Distanza euclidea media dal baricentro della base, pesata con Edep |
| `generated_scint`, `generated_cherenkov`, `injected_optical` | I | Fotoni effettivamente creati nel BGO, separati per origine; injected serve ai test |
| `arrived_scint`, `arrived_cherenkov`, `arrived_injected` | I | Primi arrivi sulla base **intera**, prima di qualsiasi selezione |
| `bulk_absorbed`, `surface_absorbed`, `escaped` | I | Perdite nel volume, sulla vernice/superficie e trasmissioni verso l'esterno |
| `reflections`, `diffuse_reflections`, `specular_reflections` | I | Conteggi aggregati dell'evento; il totale include anche TIR/Fresnel nei confronti senza vernice |
| `max_optical_steps` | I | Massimo numero di passi di una traccia ottica |
| `boundary_repairs`, `secondary_birth_repairs` | I | Contatori diagnostici delle correzioni numeriche |
| `max_boundary_shift_mm`, `max_birth_shift_mm` | D | Massimo spostamento applicato |

Con Edep=0, le due distanze pesate valgono NaN. Gli eventi restano nel file.
L'energia depositata dall'assorbimento di fotoni ottici **non** è aggiunta a Edep,
per evitare di contare due volte l'energia sorgente.

Bilancio verificato per ogni evento:

```text
generated_scint + generated_cherenkov + injected_optical
 = arrived_scint + arrived_cherenkov + arrived_injected
 + bulk_absorbed + surface_absorbed + escaped
```

Non sono salvati rapporti di analisi nella produzione. Per una selezione,
`somma(N_flux_selezionati, origin=1)/somma(generated_scint)` dà l'efficienza,
`somma(N_flux_selezionati, origin=1)/somma(edep_MeV)` dà fotoni/MeV.
Dividendo quest'ultimo per la resa impostata si ottiene la normalizzazione alla
resa di riferimento. Usare tutti gli eventi nel denominatore, inclusi quelli senza arrivi.

## `BDX::Edep` — item 3

| Campi | Tipo | Significato/unità |
|---|---|---|
| `source_step` | I | Indice della riga in questa banca, da zero per evento |
| `track`, `parent`, `pid`, `step` | I | ID traccia, genitore, codice PDG e numero di passo Geant4 della particella sorgente |
| `x0_mm`, `y0_mm`, `z0_mm`, `x1_mm`, `y1_mm`, `z1_mm` | D | Estremi iniziale/finale del passo, coordinate globali |
| `t0_ns`, `t1_ns` | D | Tempi globali degli estremi |
| `edep_MeV` | D | Deposizione non ottica nel passo |
| `generated_scint`, `generated_cherenkov` | I | Numero effettivo di fotoni creati in questo passo |

Somma Edep e conteggi dei passi coincidono con `BDX::Event`. Un passo senza Edep
ma con produzione Cherenkov è incluso. Le posizioni di nascita dei fotoni raccolti
sono salvate individualmente in Flux; non si assume che nascano tutti nel punto medio.

## `BDX::Flux` — item 4

| Campi | Tipo | Significato/unità |
|---|---|---|
| `track`, `parent` | I | ID della traccia ottica e del genitore, unici nel singolo evento |
| `origin` | I | **1 scintillazione**, **2 Cherenkov**, **0 primario ottico di test** |
| `source_step` | I | Riga Edep del passo che ha creato il fotone; -1 per primari ottici |
| `x_mm`, `y_mm`, `z_mm` | D | Primo impatto sulla base maggiore, coordinate globali; z≈120 mm |
| `lx_mm`, `ly_mm` | D | Coordinate sul piano rispetto al baricentro: lx=x, ly=y-center_y |
| `ux`, `uy`, `uz` | D | Direzione unitaria incidente nel BGO, prima della frontiera; uz>0 |
| `time_ns` | D | Tempo globale all'arrivo, include tempo di emissione e propagazione |
| `energy_eV`, `wavelength_nm` | D | Energia e lunghezza d'onda ottica |
| `path_mm` | D | Cammino totale del fotone fino al piano di raccolta |
| `steps`, `reflections` | I | Passi e riflessioni di questo fotone prima dell'arrivo |
| `birth_x_mm`, `birth_y_mm`, `birth_z_mm`, `birth_time_ns` | D | Posizione/tempo alla creazione, prima delle eventuali correzioni numeriche |
| `weight` | D | Peso traccia, sempre 1 nelle produzioni attuali |

Ogni fotone arriva **una sola volta**: subito dopo il conteggio viene fermato.
Non è applicato un taglio di area, tempo, energia, angolo o PDE.
Non sono salvati singolarmente i fotoni persi: i denominatori sono in Event/Edep.
Il collettore ha spessore 0.1 mm ma si registra l'ingresso a +120 mm, non l'uscita.

Esempi di selezioni sul CSV esportato:

```python
import csv
from collections import Counter

counts = Counter()
with open('Flux.csv') as f:
    for photon in csv.DictReader(f):
        if int(photon['origin']) != 1:
            continue
        if abs(float(photon['lx_mm'])) <= 1.5 and abs(float(photon['ly_mm'])) <= 1.5:
            counts[int(photon['record_index'])] += 1

# Leggere tutti gli eventi da Event.csv; counts[event] restituisce 0 se assente.
# Per la base intera rimuovere il solo taglio lx/ly.
```

Un taglio a posteriori non ripristina la vernice fuori da un piccolo SiPM:
è l'accettanza dei primi arrivi su una base interamente assorbente.

## Banche standard GEMC

La produzione usa anche il writer standard GEMC: `RUN::config` identifica
gli eventi e `MC::Particle` contiene le particelle primarie (`pid`, impulsi
in **GeV/c**, vertici in **cm**, `vt` in ns). La calibrazione e l’output del raster CLAS12 sono disabilitati per questo studio
isolato: `RASTER::adc` può esistere nel dizionario ma non è popolata.
Le banche BDX usano invece mm e MeV,
con unità esplicite nei nomi. `RUN::rf`, `MC::Lund`, `MC::Event`, `MC::True` e
le altre definizioni standard possono essere presenti nel dizionario anche
quando non contengono righe utili allo studio. `banks.json` riporta le righe
realmente lette per ciascuna banca; non confondere uno schema con una banca popolata.

Il collettore usa un contatore ottico dedicato e produce `BDX::Flux`, **non**
`FLUX::adc`: quella banca standard digitalizzata non conserva lo spazio delle
fasi dei singoli fotoni. `MC::True` non è la sorgente dei depositi BDX: usare Edep.
