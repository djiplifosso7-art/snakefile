# Documentation des données du projet HOLOSTRESS

## Fichiers de configuration

### config.yaml
Fichier principal de configuration du pipeline Snakemake.
Il contient les chemins vers les fichiers d’entrée, la référence génomique, les paramètres modkit et les ressources utilisées par les règles.

### samples.csv
Table décrivant les échantillons analysés. Chaque ligne correspond à un échantillon Nanopore.
Colonnes principales :
- sample : identifiant de l’échantillon
- chemin : chemin vers le fichier BAM brut ou prétraité
- temperature : température expérimentale
- stress : statut de stress thermique
- symbiont : statut symbiotique si disponible
- generation : génération si disponible

### comparaisons.csv
Table décrivant les comparaisons expérimentales utilisées pour les analyses différentielles.
Chaque ligne correspond à une comparaison entre deux groupes d’échantillons.
Colonnes principales :
- comparaison : nom de la comparaison
- samples_A : liste des échantillons du groupe A, séparés par des points-virgules (groupe A et B ce sont les différentes conditions que je compare)
- samples_B : liste des échantillons du groupe B, séparés par des points-virgules 

## Références génomiques

### aphid.fna
Génome de référence du puceron du pois Acyrthosiphon pisum.

### serratia.fna
Génome de référence de Serratia symbiotica.

### mito_aphid.only.fa
Séquence mitochondriale du puceron.

### combined_with_mito_new.fa
Référence combinée utilisée pour l’alignement des reads Nanopore.
Elle contient les séquences du puceron, de Serratia et de la mitochondrie.
###########################################################################
## Localisation des données et des fichiers intermédiaires
Afin de faciliter la reproductibilité du pipeline et d’éviter de surcharger le dossier principal du projet, les fichiers ont été organisés 
selon leur nature et leur taille. Les données brutes et les fichiers intermédiaires volumineux sont stockés dans des espaces dédiés du cluster,
tandis que les résultats légers, les tables, les figures et les logs sont conservés dans le dossier principal du
 pipeline (répertoire Pipeline_methylation).
### Données brutes issues du séquenceur Nanopore
Les fichiers issus du séquençage Nanopore sont localisés dans les dossiers propres à chaque échantillon :
```bash
/env/cns/proj/projet_DWS/{sample}/RunsNanopore/
#Dans le pipeline, les fichiers BAM sont recherchés automatiquement dans ce dossier pour chaque échantillon. Les fichiers contenant 
#_fail ainsi que les dossiers CLONE sont exclus afin de ne conserver que les fichiers utilisables pour l’analyse.
### Dossiers important pour le pipeline
##Le pipeline est localisé dans : /env/cns/proj/projet_DWS/Pipeline_methylation/
Ce dossier contient les fichiers de configuration, les métadonnées et les résultats finaux légers :
config.yaml
samples.csv
comparaisons.csv
snakefile
results/
logs/
##Fichiers intermédiaires volumineux
Les fichiers intermédiaires les plus lourds sont stockés dans l’espace temporaire BIGTMP :
/env/cns/bigtmp/adjiplif/HOLOSTRESS/results/{sample}/{sample}.aligned.sorted.bam
/env/cns/bigtmp/adjiplif/HOLOSTRESS/results/{sample}/{sample}.aligned.sorted.bam.bai
/env/cns/bigtmp/adjiplif/HOLOSTRESS/results/{sample}/{sample}.{organism}.aligned.sorted.bam
/env/cns/bigtmp/adjiplif/HOLOSTRESS/results/{sample}/{sample}.{organism}.aligned.sorted.bam.bai
/env/cns/bigtmp/adjiplif/HOLOSTRESS/results/{sample}/{sample}.{organism}.{modtype}.mods.raw.bed
/env/cns/bigtmp/adjiplif/HOLOSTRESS/results/{sample}/{sample}.{organism}.{modtype}.mods.bed

Ces fichiers sont volumineux car ils contiennent les reads alignés ou les informations de méthylation à l’échelle du génome. Ils sont nécessaires 
au fonctionnement du pipeline, mais ne correspondent pas directement aux résultats finaux à interpréter.
##Fichiers intermédiaires compressés et indexés
Les fichiers bedMethyl filtrés sont également compressés et indexés afin d’être utilisés par certains outils, notamment modkit dmr pair
. Ils sont stockés dans l’espace SCRATCH :
/env/cns/proj/projet_DWS/scratch/adjiplif/HOLOSTRESS/results/

##Fichiers intermédiaires légers et résultats finaux
Les fichiers moins volumineux sont conservés dans le dossier principal du pipeline :
/env/cns/proj/projet_DWS/Pipeline_methylation/results/
Ils comprennent notamment les contrôles qualité, les tables de résultats, les figures, les fichiers convertis pour methylKit,
 les annotations et les logs.
results/qc/
results/methylkit/input/
results/methylkit/diff/
results/methylkit/tables/
results/methylkit/plots/
results/methylkit/annotation/
results/methylkit/enrichment/
results/methylkit/logs/

##Fichiers methylKit
Les fichiers d’entrée pour methylKit sont stockés dans :
results/methylkit/input/ 
exple: results/methylkit/input/{sample}.{organism}.{modtype}.methylkit.tsv
Les résultats de méthylation différentielle sont stockés dans :
results/methylkit/diff/
Les tables et figures finales sont stockées dans :
results/methylkit/tables/
results/methylkit/plots/
Les annotations des DMR et les enrichissements fonctionnels sont stockés dans :
results/methylkit/annotation/
results/methylkit/enrichment/

#####  SORTIES IMPORTANTES
# Organisation des résultats

## results/qc/
Contient les figures et tables de contrôle qualité :
- qualité avant alignement
- taux d’alignement
- couverture par organisme
- qualité avant détection des DMR

## results/methylkit/input/
Contient les fichiers convertis au format compatible avec methylKit.
Ces fichiers sont produits à partir des fichiers bedMethyl générés par modkit.

Format attendu :
chr | start | end | strand | coverage | freqC

## results/methylkit/diff/
Contient les résultats des comparaisons différentielles methylKit :
- *.counts.tsv : résumé du nombre de DMR hyper- et hypométhylées
- *.hyper.tsv : DMR hyperméthylées
- *.hypo.tsv : DMR hypométhylées

## results/methylkit/tables/
Contient les tables récapitulatives utilisées pour les figures.

## results/methylkit/plots/
Contient les figures finales issues de methylKit, notamment les histogrammes du nombre de DMR hyper- et hypométhylées.

## results/methylkit/annotation/
Contient les DMR annotées avec les gènes proches ou chevauchants.

## results/methylkit/enrichment/
Contient les résultats d’enrichissement fonctionnel GO réalisés à partir des gènes associés aux DMR.

## logs/
Contient les logs généraux du pipeline et des jobs SLURM.
+999

