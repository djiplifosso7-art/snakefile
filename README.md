# Pipeline d’analyse du méthylome Nanopore - Projet HOLOSTRESS

Ce pipeline Snakemake permet d’analyser les données de méthylation issues du séquençage Oxford Nanopore dans le cadre du projet HOLOSTRESS.

L’objectif est d’identifier des différences de méthylation associées au stress thermique, à la température d’acclimatation,
 à la présence de symbiotes et à la génération chez le puceron du pois Acyrthosiphon pisum.

## Étapes principales du pipeline

1. Fusion des fichiers BAM bruts par échantillon
2. Alignement des reads Nanopore sur une référence combinée
3. Tri et indexation des fichiers BAM
4. Séparation des BAM par organisme
5. Appel des bases modifiées avec modkit
6. Filtrage des fichiers bedMethyl selon la couverture
7. Conversion des fichiers modkit vers un format compatible methylKit
8. Détection des régions différentiellement méthylées avec methylKit
9. Génération des figures de synthèse
10. Annotation des DMR aux gènes
11. Enrichissement fonctionnel GO

## Fichiers de configuration

Le pipeline utilise trois fichiers principaux :

- `config.yaml` : paramètres généraux du pipeline
- `samples.csv` : description des échantillons
- `comparaisons.csv` : description des comparaisons expérimentales

## Lancer un dry-run

```bash
snakemake -s snakefile -np

