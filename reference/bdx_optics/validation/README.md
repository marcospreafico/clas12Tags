# Verifiche numeriche del modello ottico BDX

- 6 controlli Python, incluso confronto di tutti i vertici con il Perl originale.
- 15 controlli Geant4: raggi diretti, apertura del sensore, riflettivita, assorbimento,
  riflessione totale interna nel caso senza vernice, tre zone della vernice,
  raggi oltre 100 passi e convergenza con boundary_push_mm = 1e-7, 1e-6, 1e-5.
- Scansione degli otto tipi: 10, 120 e 230 mm dalla faccia SiPM; 3 muoni/punto;
  resa ridotta a 100 fotoni/MeV; 72 eventi, 277651 fotoni generati.
- Controllo alla resa completa: tipo 8, 120 mm, 2 muoni, 8500 fotoni/MeV;
  658585 fotoni generati e 6202 arrivati. Resa generata misurata: 8511.84 fotoni/MeV;
  raccolta = 0.00941716; fotoni arrivati/MeV = 80.1574.

Tutti gli eventi chiudono esattamente il bilancio dei fotoni. I controlli delle
probabilita sono confronti statistici, non uguaglianze dei singoli conteggi
alle loro medie teoriche. Gli esiti e i conteggi delle prime riflessioni sono
in transport_checks.json; le configurazioni in transport_inputs.json.

Le curve eight_types_collection.png/svg sono una verifica preliminare con poca
statistica e solo tre distanze, non risultati di produzione. Le linee collegano
semplicemente i punti e non rappresentano un fit. Per confronti precisi aumentare
statistica e punti lungo z usando scripts/run_bdx_optics.py.

I manifest e i CSV degli eventi sono conservati qui; i run completi con input,
log e arrivi individuali restano in build/bdx_optics/eight_types_check e
build/bdx_optics/full_yield_check. Le due campagne registrano i sorgenti e i semi
usati alla partenza; le analisi registrano anche l'hash della propria versione.
I conteggi di validazione non validano sperimentalmente il modello della vernice.

Per la correzione successiva dei casi limite di navigazione e i risultati aggiornati, vedere [la regressione del 24 novembre 2026](20261124_boundary_fix/README.md).
