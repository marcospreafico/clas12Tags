# Regressione del run fallito: 24 novembre 2026

La correzione conserva vertici, normali e superfici del GenericTrap originale.
Rende conservative le stime scalari della distanza dalla parete e risolve le
intersezioni dei raggi con formule quadratiche stabili. La piccola correzione dei punti di
riflessione/nascita rimane limitata e verificata; la generazione di scintillazione
rimane nativa, senza ricampionamento dei punti di emissione.

Verifiche completate:

- Il punto originale BGO-1, distanza 20 mm, seme 2010314138: 10 eventi, resa 800 fotoni/MeV.
- Tutte le 32 combinazioni del comando utente: otto tipi e distanze 20/80/140/200 mm,
  con un evento per punto e la stessa resa di 800 fotoni/MeV.
- 48067 punti sulle otto geometrie: caso numerico originale, superfici, spigoli e
  punti esterni; rifiuto di correzioni macroscopiche e verifiche della distanza conservativa.
- 16 controlli Geant4 di trasporto, compreso il test di fotoni generati al bordo;
  sei controlli Python di geometria e analisi.

Tutti i bilanci dei fotoni dei run di verifica chiudono esattamente. La verifica
sulle 32 combinazioni usa un evento per punto; non costituisce una scansione di
produzione con dieci eventi per punto. I conteggi e gli hash sono in summary.json.
I run originali falliti restano separati nei rispettivi percorsi in build/bdx_optics.
