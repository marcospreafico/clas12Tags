# Provenienza e ipotesi: BDX BGO

## Geometria documentata

- Repository: <https://github.com/marcospreafico/clas12Tags>, ramo `frederic`.
- Commit clonato: `1d2313f8cc701e43bf461bee44ea12f43c7d402d`.
- File: `experiments/bdx/geometry_crs60.pl`.
- SHA256: `04d3f3284707593209545b2a9a41a54326736258ff212ee2d571dd0b34ba0218`.
- Otto definizioni attive di `G4GenericTrap`; conservati tutti i vertici e il
  loro ordine. `theta` e presente nella tabella ma non usato dalla costruzione
  attiva originale: non viene introdotta una rotazione aggiuntiva.
- Le dimensioni Perl sono semidimensioni in cm; estrazione in mm tramite
  `scripts/bdx_geometry.py`. Lunghezza totale 240 mm per ogni tipo.
- Il test esegue il Perl originale con una funzione di cattura `print_det` e
  confronta numericamente le 17 dimensioni emesse per ognuno degli otto tipi.

Derivazione del centro del SiPM: per la base maggiore trapezoidale, simmetrica
rispetto a x, con semibasi `Bm` a y=-H e `Bg` a y=+H, il baricentro di area e
`x=0`, `y=H*(Bg-Bm)/(3*(Bg+Bm))`. Le coordinate locali del flux detector e le selezioni di analisi usano questo centro.
Il caso tipo 8, con basi uguali, ha y=0.

## Dati e scelte dell'utente (23 novembre 2026)

- Richiesta iniziale: sensore 3 x 3 mm² al centro della base maggiore lungo z.
- Aggiornamento del 24 novembre 2026: flux detector su tutta la base maggiore,
  primo arrivo salvato per fotone, selezioni spaziali solo in analisi.
- Lunghezza di assorbimento BGO: **10 m = 10000 mm**, comunicata il 24 novembre 2026.
- Accoppiamento ottico da trascurare: ingresso ideale con indice adattato al BGO.
- Rivestimento corretto: **vernice bianca**, non mylar.
- R assoluta = 0.98 a ogni interazione.
- Da base maggiore: 40 mm diffusi opachi, 100 mm di miscela continua,
  ultimi 100 mm speculari. Tutte le lunghezze sono configurabili.

Interpretazione modellistica: "diffusore opaco" = Lambertiano riflettente,
"riflettente" = speculare. La frazione speculare aumenta linearmente da 0 a 1
nella zona di transizione. Le probabilita non condizionate sono
`P_ass=1-R`, `P_spec=R*f_spec(d)`, `P_diff=R*(1-f_spec(d))`.
La base minore e le facce laterali seguono il modello alla loro z. La base
maggiore intera è un collettore ideale aperto; la vernice non vi è applicata.

Implementazione: modello Geant4 UNIFIED `dielectric_metal`, finitura `ground`,
`sigma_alpha=0`, `SPECULARLOBECONSTANT=BACKSCATTERCONSTANT=0`,
`SPECULARSPIKECONSTANT=f_spec(d)` prima del campionamento nativo della frontiera.
Il complemento produce riflessioni Lambertiane. R e `REFLECTIVITY`; `EFFICIENCY=0`
impedisce di confondere una perdita di vernice con una rivelazione nel SiPM.
Il volume del cristallo resta un solo solido continuo.

I riferimenti sul mylar consultati prima della correzione dell'utente non sono
utilizzati nei parametri attivi o nei risultati della simulazione.

## Riferimento BGO indipendente da MINE

[Hilger Crystals, Bismuth Germanate](https://www.hilger-crystals.co.uk/materials/bismuth-germanate/)
riporta 8500 fotoni/MeV, picco 480 nm, indice al picco 2.15 e decadimento 300 ns.
Questi valori sono un riferimento a temperatura ambiente, non misure dei
cristalli BDX. La composizione e la densita per il trasporto della particella
sono quelle del materiale NIST `G4_BGO` della versione Geant4 registrata nel log.

Ipotesi per isolare la geometria: indice costante, emissione a picco stretto,
assorbimento nel volume con lunghezza costante 10000 mm, scintillazione proporzionale senza Birks,
nessuna diffusione Rayleigh/Mie/WLS, Cherenkov disattivato. Eventuale Cherenkov
abilitato dall'utente resta separato dalla normalizzazione di scintillazione.
L'emissione stretta non ricostruisce lo spettro misurato del BGO; con coefficienti
indipendenti dalla lunghezza d'onda non modifica la raccolta ottica integrata.

Guida primaria per i processi ottici:
[Geant4 Physics Processes](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/TrackingAndPhysics/physicsProcess.html).
Il codice gestisce i nomi delle proprieta di scintillazione Geant4 10.7 e 11;
la compilazione e le verifiche locali sono eseguite con **10.7.4**. La versione
11 non e stata verificata in questo ambiente.

## Output e verifica numerica

`N_generated` conta i secondari effettivi creati da `Scintillation`; `N_arrived`
conta solo ingressi realmente trasmessi nel collettore sulla base intera. I fotoni riflessi non vengono
scambiati per ingressi, anche se il volume candidato della navigazione e il collettore.
Ogni traccia arrivata viene fermata una volta; la PDE non e applicata.

Bilancio obbligatorio evento per evento:

```text
generati_scint + generati_Cherenkov + primari_ottici
 = arrivati_scint + arrivati_Cherenkov + arrivati_primari
 + assorbiti_volume + assorbiti_superficie + usciti
```

Nessun taglio su energia, tempo o numero di riflessioni dei fotoni. Geant4 10.7.4
richiede una correzione numerica di posizione/touchable dopo le riflessioni sul
`G4GenericTrap`: si corregge l'eventuale distanza di arrotondamento fuori dalla
superficie e si aggiunge una distanza interna di riferimento di 1 nm; vicino a
spigoli/facce non piane si verifica un punto interno cercato verso il centro.
La correzione e limitata a 1 micrometro piu la distanza interna configurata.
Si sincronizzano volume e velocita. Non si cambia direzione, energia o probabilita di riflessione.
Questa correzione viene verificata con raggi analitici e variazioni del parametro
`boundary_push_mm`. Un transito finito attribuito all'aria all'interno del BGO
fa fallire il run.

Le verifiche numeriche confermano l'implementazione; non sono una validazione
sperimentale di R, della miscela angolare o della resa dei cristalli BDX.


## Correzione del caso limite del 24 novembre 2026

Queste regressioni storiche precedono il collettore sulla base intera e usavano
SiPM 3×3 mm², assorbimento nullo. I risultati conservati non vanno reinterpretati
come nuove produzioni con assorbimento di 10 m.


Il run BGO-1, distanza 20 mm, resa 800 fotoni/MeV e seme 2010314138 ha prodotto
un'intersezione a (-14.554468366128804, -22.489342449694302, 75.987450886933132) mm.
Geant4 classifica questo punto come esterno per circa 1.276 nm. Il precedente
spostamento fisso di 1 nm era insufficiente e interrompeva l'evento. Il nuovo
riposizionamento verifica la distanza numerica residua e l'appartenenza al solido,
con un limite esplicito che impedisce di correggere errori macroscopici.

`tests/bdx_optics/relocation_test.cc` contiene quel punto di regressione,
campiona tutte le otto geometrie (superfici, punti spostati fuori e vertici) e
controlla che punti macroscopicamente esterni siano rifiutati. I CSV di evento
registrano numero di riparazioni e spostamento massimo; l'analisi li aggrega.


La ripetizione estesa ha isolato un secondo problema, nella distanza minima
usata dalla navigazione e dalla diffusione multipla su facce non piane. Nel
punto (-20.994472375460159,-0.027083894775988219,100.03453154064427) mm del tipo 1,
`G4GenericTrap::DistanceToOut(p)` restituisce 0.0257894 mm, mentre esiste un punto
sulla parete a distanza inferiore a 0.00683104 mm. La sovrastima puo far uscire
un elettrone dal BGO mantenendo il materiale del passo. Il problema e riprodotto
in `relocation_test.cc` e confrontato con il codice di
[Geant4 10.7.4 G4GenericTrap](https://github.com/Geant4/geant4/blob/v10.7.4/source/geometry/solids/specific/src/G4GenericTrap.cc).

La classe locale `ConservativeBgoTrap` eredita il solido originale, conservando
vertici, superfici, normali e classificazione di appartenenza. Sostituisce le
distanze scalari con limiti inferiori conservativi e risolve le intersezioni
direzionali sulle stesse superfici con formule numericamente stabili.
Per ogni z si scrive la sezione come `-h(z)<=y<=h(z)` e
`|x|<=a(z)+c(z)*y/h(z)`, dove h, a e c interpolano linearmente i valori dei
vertici. Le due superfici laterali sono descritte da
`F=h(z)*(a(z) +/- x)+c(z)*y`. Il rapporto fra il margine F e un limite superiore
della norma del suo gradiente nel bounding box fornisce un limite inferiore alla
distanza dalla superficie. Per le facce y/z si usano i corrispondenti margini
normalizzati; fuori dal bounding box la distanza da esso e gia un limite valido.
Nel punto catturato la nuova stima e 0.00644736 mm, inferiore al limite geometrico.

I test confrontano queste stime con distanze da punti noti delle superfici di
tutte le otto forme e verificano che la geometria rimanga identica. La generazione
Geant4 dei fotoni e mantenuta senza ricampionamento delle posizioni: viene solo
riparato l'arrotondamento dei fotoni nati sul bordo, prima del trasporto.
Il test `secondary_birth_outside_BGO` forza questa condizione in dieci eventi
senza eliminare o ripesare fotoni e verifica il bilancio completo.


Le intersezioni usano gli stessi sei vincoli del volume: quattro lineari e due
quadratici lungo un raggio. Le radici quadratiche sono calcolate con la forma
stabile `q=-0.5*(b+copysign(sqrt(discriminante),b))`, `q/a` e `c/q`, evitando la
cancellazione fra termini quasi uguali. Il segno della derivata determina i
passaggi in ingresso e in uscita; tutti i vincoli devono essere soddisfatti al
punto di attraversamento. I test comprendono raggi trasversali a quattro z
per ciascun tipo, componenti z nulle/quasi nulle e raggi casuali con controllo
di uscita e rientro inverso. Questo intercetta anche attraversamenti falsamente
quasi nulli, che produrrebbero energia depositata o fotoni praticamente zero.
