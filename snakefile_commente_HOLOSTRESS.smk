# =============================================================================
# Snakefile - Pipeline d'analyse de méthylation ONT pour le projet HOLOSTRESS
# =============================================================================
# Objectif général :
#   1) Lire la configuration du projet et les métadonnées des échantillons.
#   2) Fusionner, aligner, trier et indexer les BAM Nanopore.
#   3) Séparer les alignements par organisme à partir d'une référence combinée.
#   4) Appeler les marques de méthylation avec modkit : 4mC, 5mC et 6mA.
#   5) Convertir les sorties modkit vers un format compatible methylKit.
#   6) Détecter des régions/sites différentiellement méthylés avec methylKit.
#   7) Résumer les DMR hyper/hypométhylées, annoter les DMR aux gènes
#      et réaliser l'enrichissement fonctionnel.
#
# Remarque :
#   Les commentaires ajoutés ici visent à rendre le workflow lisible pour une
#   personne qui n'a pas écrit le code. Ils n'ont pas vocation à modifier la
#   logique du pipeline.
# =============================================================================

# Fichier de configuration principal : contient les chemins vers les données,
# les fichiers CSV de métadonnées, la référence et les paramètres des outils.
configfile: "config.yaml"
# pandas est utilisé ici pour lire les tables CSV de métadonnées.
import pandas as pd

# Table des échantillons : doit contenir au moins la colonne "sample"
# et les chemins vers les BAM bruts.
samplesheet = pd.read_csv(config["samples"])
# Table des comparaisons biologiques : chaque ligne définit un contraste,
# avec les échantillons du groupe A et ceux du groupe B.
comparaisons = pd.read_csv(config["comparaisons"])

# Liste des échantillons à traiter dans le workflow.
SAMPLES = samplesheet["sample"].tolist()
# Liste des comparaisons demandées pour les analyses différentielles.
COMPARAISONS = comparaisons["comparaison"].tolist()
comparison_map = comparaisons.set_index("comparaison").to_dict("index")
# Génome de référence combiné utilisé pour l’alignement.
# Il contient le puceron, Serratia et la mitochondrie du puceron.
REF = config["reference"]

# Organismes analysés séparément après alignement sur la référence combinée.
ORGANISMS = ["aphid", "serratia"]

BIGTMP = "/env/cns/bigtmp/adjiplif/HOLOSTRESS"
SCRATCH = "/env/cns/proj/projet_DWS/scratch/adjiplif/HOLOSTRESS"
BIGTMP_RESULTS = f"{BIGTMP}/results"
SCRATCH_RESULTS = f"{SCRATCH}/results"

CONTIG_LISTS = {
    "aphid": "contig_aphid_list",
    "serratia": "contig_serratia_list",
}


# Retourne la référence FASTA spécifique à l’organisme demandé.
# Cette fonction est utilisée par les règles qui doivent travailler organisme par organisme.
def organism_ref(wildcards):
    return f"{SCRATCH}/refs/{wildcards.organism}.only.fa"
	
# Marques de méthylation étudiées.
MOD_TYPES = ["4mC", "5mC", "6mA"]


MOD_BASES = {
    "4mC": "C:21839",
    "5mC": "C:m",
    "6mA": "A:a",
}

DMR_BASE = {
    "4mC": "C",
    "5mC": "C",
    "6mA": "A",
}

# Récupère dans la table des échantillons le chemin du BAM brut associé à un sample.
def get_bam_path(wildcards):
    row = samplesheet[samplesheet["sample"] == wildcards.sample]
    return row.iloc[0]["chemin"]

MODKIT_EXTRA = config["modkit"].get("extra", "")
MODKIT_MIN_COV = config["modkit"].get("min_coverage", 10)
##Verifier que les samples du csv correspondenty à ceux de mon config
#config_samples = set(config["samples"].keys())
#csv_samples = set(samplesheet["sample"].tolist())

#missing_in_config = csv_samples - config_samples
#missing_in_csv = config_samples - csv_samples

#assert not missing_in_config, f"Samples présents dans samples.csv mais absents du config.yaml: {missing_in_config}"
#assert not missing_in_csv, f"Samples présents dans config.yaml mais absents du samples.csv: {missing_in_csv}"

##récupérer un métadata depuis le csv
# Fonction générique pour extraire une information de métadonnées
# depuis le fichier samples.csv.
def get_sample_info(wildcards, column):
    row = samplesheet[samplesheet["sample"] == wildcards.sample]
    return row.iloc[0][column]

def get_temperature(wildcards):
    return get_sample_info(wildcards, "temperature")

def get_stress(wildcards):
    return get_sample_info(wildcards, "stress")
###Mes fonctions pour le dmr pair de mes conditions
##fonction pour extraire les samples d'une comparaison
# Pour une comparaison donnée, récupère les échantillons du groupe A
# et ceux du groupe B depuis comparaisons.csv.
def get_samples_for_comparison(comp):
    row = comparison_map[comp]
    group_a = [x.strip() for x in str(row["samples_A"]).split(";") if x.strip()]
    group_b = [x.strip() for x in str(row["samples_B"]).split(";") if x.strip()]

  #  print(f"[DEBUG] {comp}: group_a={group_a}, group_b={group_b}")
    return group_a, group_b

## Ajouter une liste de comparaisons valides
# Vérifie qu’une comparaison est utilisable :
#   - les deux groupes contiennent au moins un échantillon ;
#   - tous les échantillons cités existent dans samples.csv.
def comparison_is_valid(comp):
    group_a, group_b = get_samples_for_comparison(comp)
    return (
        len(group_a) > 0 and len(group_b) > 0
        and all(s in SAMPLES for s in group_a)
        and all(s in SAMPLES for s in group_b)
    )

VALID_COMPARAISONS = [
    comp for comp in COMPARAISONS
    if comparison_is_valid(comp)
]

DMR_ORG_MODS = [("aphid","4mC"), ("aphid","5mC"), ("aphid","6mA")]

DMR_TARGETS = [
    f"results/dmr_pair/{comp}.{org}.{mod}.bed"
    for comp in VALID_COMPARAISONS
    for org, mod in DMR_ORG_MODS
]


# ============================================================
# methylKit downstream analysis
# ============================================================

# Dossier central pour les fichiers d’entrée, sorties différentielles,
# figures et annotations liés à methylKit.
METHYLKIT_DIR = "results/methylkit"

# On teste toutes les marques et organismes disponibles
METHYLKIT_ORG_MODS = [
    (org, mod)
    for org in ORGANISMS
    for mod in MOD_TYPES
]

METHYLKIT_INPUTS = [
    f"{METHYLKIT_DIR}/input/{sample}.{org}.{mod}.methylkit.tsv"
    for sample in SAMPLES
    for org, mod in METHYLKIT_ORG_MODS
]

METHYLKIT_COUNTS = [
    f"{METHYLKIT_DIR}/diff/{comp}.{org}.{mod}.counts.tsv"
    for comp in VALID_COMPARAISONS
    for org, mod in METHYLKIT_ORG_MODS
]

# ============================================================
# Annotation et enrichissement des DMR methylKit
# ============================================================

GENE_BEDS = {
    "aphid": "results/annotation/aphid.genes.bed",
    "serratia": "results/annotation/serratia.genes.bed",
}

GENE2GO = {
    "aphid": "results/annotation/aphid_gene2go.tsv",
    # "serratia": "results/annotation/serratia_gene2go.tsv",
}

# Distance autour des gènes pour associer une DMR à un gène.
# 0 = uniquement chevauchement direct avec un gène.
# 2000 = gènes chevauchants ou situés à moins de 2 kb.
DMR_GENE_WINDOW = 2000

METHYLKIT_DMR_ANNOTATIONS = [
    f"{METHYLKIT_DIR}/annotation/{comp}.{org}.{mod}.dmr_annotated.tsv"
    for comp in VALID_COMPARAISONS
    for org, mod in METHYLKIT_ORG_MODS
]

METHYLKIT_DMR_GENE_LISTS = [
    f"{METHYLKIT_DIR}/annotation/{comp}.{org}.{mod}.dmr_genes.txt"
    for comp in VALID_COMPARAISONS
    for org, mod in METHYLKIT_ORG_MODS
]

# Pour l’enrichissement, je te conseille de commencer par aphid uniquement,
# car les DMR Serratia ne sont pas assez robustes à cause de la couverture.
METHYLKIT_GO_ENRICHMENTS = [
    f"{METHYLKIT_DIR}/enrichment/{comp}.aphid.{mod}.GO_enrichment.tsv"
    for comp in VALID_COMPARAISONS
    for mod in MOD_TYPES
]

# Construit la liste des fichiers methylKit correspondant au groupe A
# pour une comparaison, un organisme et une marque donnés.
def methylkit_a_files(wildcards):
    group_a, _ = get_samples_for_comparison(wildcards.comp)
    return [
        f"{METHYLKIT_DIR}/input/{s}.{wildcards.organism}.{wildcards.modtype}.methylkit.tsv"
        for s in group_a
    ]

# Construit la liste des fichiers methylKit correspondant au groupe B
# pour une comparaison, un organisme et une marque donnés.
def methylkit_b_files(wildcards):
    _, group_b = get_samples_for_comparison(wildcards.comp)
    return [
        f"{METHYLKIT_DIR}/input/{s}.{wildcards.organism}.{wildcards.modtype}.methylkit.tsv"
        for s in group_b
    ]

### pour créer mes input de modkit
# =============================================================================
# Règle finale du workflow
# =============================================================================
# Cette règle liste les fichiers finaux que Snakemake doit produire.
# Elle ne contient pas de commande : elle sert à définir les objectifs du pipeline.
# Si un fichier listé ici est absent, Snakemake cherche automatiquement la règle
# capable de le produire.
rule all:
    input:
        expand(f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.aligned.sorted.bam", sample=SAMPLES),
        expand(f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.aligned.sorted.bam.bai", sample=SAMPLES),
        expand(f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.flagstat.txt", sample=SAMPLES),
        expand(f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam", sample=SAMPLES, organism=ORGANISMS),
        expand(f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam.bai", sample=SAMPLES, organism=ORGANISMS),
        expand("results/{sample}/{sample}.{organism}.coverage.tsv", sample=SAMPLES, organism=ORGANISMS),
        expand("results/{sample}/{sample}.{organism}.idxstats.txt", sample=SAMPLES, organism=ORGANISMS),
        #expand("results/{sample}/{sample}.{organism}.modkit_summary.txt", sample=SAMPLES, organism=ORGANISMS),
        #expand(f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed", sample=SAMPLES, organism=ORGANISMS, modtype=MOD_TYPES),
        #expand(f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed.gz", sample=SAMPLES, organism=ORGANISMS, modtype=MOD_TYPES),
        #expand(f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed.gz.tbi", sample=SAMPLES, organism=ORGANISMS, modtype=MOD_TYPES),
        #DMR_TARGETS,
        "contig_aphid_list",
        "contig_serratia_list",
        "contig_mito_aphid_list",
        METHYLKIT_INPUTS,
        METHYLKIT_COUNTS,
        f"{METHYLKIT_DIR}/tables/methylkit_hyper_hypo_counts_all.tsv",
        f"{METHYLKIT_DIR}/plots/methylkit_hyper_hypo_all_comparisons.png",
        f"{METHYLKIT_DIR}/plots/methylkit_hyper_hypo_all_comparisons_scale_6mA.png",
        METHYLKIT_DMR_ANNOTATIONS,
        METHYLKIT_DMR_GENE_LISTS,
        #METHYLKIT_GO_ENRICHMENTS,
        #"results/qc/qc_before_align.png",
       #"results/qc/qc_before_pileup.png",
        #"results/qc/qc_coverage_by_organism.png",
        #"results/qc/qc_before_dmr.png",
        #"results/global_methylation/global_methylation_from_beds.tsv",
        #"results/global_methylation/global_methylation_summary.tsv",
        #"results/global_methylation/global_methylation_violin.png"       



# =============================================================================
# Règle groupée pour produire les figures de contrôle qualité
# =============================================================================
# Cette règle sert de raccourci pour générer les principales figures QC.
rule qc:
    input:
        "results/qc/qc_before_align.png",
        "results/qc/qc_before_pileup.png",
        "results/qc/qc_before_dmr.png",
        "results/qc/qc_coverage_by_organism.png"



# Recherche automatiquement tous les BAM bruts associés à un échantillon.
# Les fichiers issus de dossiers CLONE ou marqués fail sont ignorés.
def raw_bams_for_sample(wildcards):
    import os

    base = f"/env/cns/proj/projet_DWS/{wildcards.sample}/RunsNanopore"

    bams = []

    for root, dirs, files in os.walk(base, followlinks=True):
        if "/CLONE" in root:
            continue

        for fn in files:
            if not fn.endswith(".bam"):
                continue

            if "_fail" in fn:
                continue

            if f"DWS_{wildcards.sample}_" not in fn:
                continue

            bams.append(os.path.join(root, fn))

    bams = sorted(bams)

    if len(bams) == 0:
        raise ValueError(f"Aucun BAM brut trouvé pour {wildcards.sample} dans {base}")

    return bams


# =============================================================================
# Fusion des BAM bruts d’un même échantillon
# =============================================================================
# Certains échantillons peuvent avoir plusieurs fichiers BAM issus de plusieurs
# runs Nanopore. Cette règle les regroupe en un seul BAM brut par échantillon.
# Si un seul BAM est trouvé, il est simplement copié.
rule bam_fusion:
    input:
        bams=raw_bams_for_sample
    output:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.merged.raw.bam"
    log:
        "results/{sample}/{sample}.bam_fusion.log"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=4 * 60
    params:
        n_bams=lambda wildcards, input: len(input.bams),
        first_bam=lambda wildcards, input: input.bams[0],
        bams_q=lambda wildcards, input: " ".join(f'"{b}"' for b in input.bams)
    shell:
        r"""
        set -euo pipefail
        module load samtools

        mkdir -p $(dirname {output.bam}) results/{wildcards.sample}

        echo "Sample: {wildcards.sample}" > {log}
        echo "Nombre de BAM trouvés: {params.n_bams}" >> {log}
        echo "BAM utilisés:" >> {log}
        printf '%s\n' {params.bams_q} >> {log}

        tmp="{output.bam}.tmp"
        rm -f "$tmp" {output.bam}

        if [ "{params.n_bams}" -eq 1 ]; then
            echo "Un seul BAM trouvé : copie simple" >> {log}
            cp "{params.first_bam}" "$tmp"
        else
            echo "Concaténation avec samtools cat" >> {log}
            samtools cat \
              -o "$tmp" \
              {params.bams_q} \
              >> {log} 2>&1
        fi

        echo "Contrôle simple du fichier produit" >> {log}

        if [ ! -s "$tmp" ]; then
            echo "ERREUR: le fichier temporaire est vide ou absent" >> {log}
            exit 1
        fi

        echo "Taille du BAM temporaire:" >> {log}
        ls -lh "$tmp" >> {log}

        mv "$tmp" {output.bam}

        echo "BAM fusionné final:" >> {log}
        ls -lh {output.bam} >> {log}
        """

		
## QC avant alignement, on utilise samtools stat
# =============================================================================
# QC avant alignement
# =============================================================================
# Cette règle mesure, à partir des BAM bruts fusionnés :
#   - la longueur moyenne des reads ;
#   - le nombre total de reads.
# Elle produit une figure de contrôle qualité avant l’alignement.
rule plot_qc_before_align:
    input:
        expand(f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.merged.raw.bam", sample=SAMPLES)
    output:
        "results/qc/qc_before_align.png"
    log:
        "results/qc/qc_before_align.log"
    run:
        import subprocess
        import re
        import pandas as pd
        import matplotlib.pyplot as plt
        import os

        os.makedirs("results/qc", exist_ok=True)

        data = []

        with open(log[0], "w") as lg:
            lg.write("QC before align on merged raw BAMs\n\n")

            for bam in input:
                sample = os.path.basename(bam).replace(".merged.raw.bam", "")

                lg.write(f"Sample: {sample}\n")
                lg.write(f"BAM: {bam}\n")

                try:
                    stats = subprocess.check_output(
                        f"module load samtools && samtools stats {bam}",
                        shell=True,
                        executable="/bin/bash",
                        stderr=subprocess.STDOUT
                    ).decode()

                    m_len = re.search(r"average length:\s+([\d\.]+)", stats)
                    mean_len = float(m_len.group(1)) if m_len else 0

                    m_reads = re.search(r"raw total sequences:\s+(\d+)", stats)
                    n_reads = int(m_reads.group(1)) if m_reads else 0

                    data.append([sample, mean_len, n_reads])

                    lg.write(f"mean_len={mean_len}\n")
                    lg.write(f"n_reads={n_reads}\n\n")

                except subprocess.CalledProcessError as e:
                    lg.write("ERROR samtools stats failed\n")
                    lg.write(e.output.decode(errors="replace"))
                    lg.write("\n\n")

                    data.append([sample, 0, 0])

        df = pd.DataFrame(data, columns=["sample", "mean_len", "n_reads"])

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].bar(df["sample"], df["mean_len"])
        axes[0].set_title("Longueur moyenne des reads")

        axes[1].bar(df["sample"], df["n_reads"])
        axes[1].set_title("Nombre de reads")

        for ax in axes:
            ax.tick_params(axis="x", rotation=90)

        plt.tight_layout()
        plt.savefig(output[0], dpi=300)	
		
# =============================================================================
# Construction de la référence combinée
# =============================================================================
# Cette règle concatène les séquences FASTA du puceron, de Serratia et de la
# mitochondrie du puceron, puis indexe la référence avec samtools faidx.
rule make_combined_reference_new:
    input:
        aphid="aphid.fna",
        serratia="serratia.fna",
        mito="mito_aphid.only.fa"
    output:
        fa="combined_with_mito_new.fa",
        fai="combined_with_mito_new.fa.fai"
    log:
        "logs/make_combined_reference_new.log"
    shell:
        r"""
        set -euo pipefail
        module load samtools

        mkdir -p logs

        echo "Création de la référence combinée" > {log}
        echo "Input aphid: {input.aphid}" >> {log}
        echo "Input serratia: {input.serratia}" >> {log}
        echo "Input mito: {input.mito}" >> {log}

        # Nettoyer les fins de lignes Windows éventuelles
        # et concaténer proprement les trois FASTA
        (
            sed 's/\r$//' {input.aphid}
            echo
            sed 's/\r$//' {input.serratia}
            echo
            sed 's/\r$//' {input.mito}
        ) > {output.fa}

        # Indexer la référence
        samtools faidx {output.fa} >> {log} 2>&1

        echo "Référence finale:" >> {log}
        ls -lh {output.fa} {output.fai} >> {log}

        echo "Nombre de contigs:" >> {log}
        grep -c '^>' {output.fa} >> {log}

        echo "Premiers contigs:" >> {log}
        grep '^>' {output.fa} | head >> {log}
        """	

#####sorted rule						
# =============================================================================
# Alignement, tri et indexation des reads Nanopore
# =============================================================================
# Cette règle aligne les BAM bruts sur la référence combinée avec Dorado,
# trie le BAM aligné avec samtools sort, puis crée l’index BAM.
# Elle produit aussi un flagstat pour contrôler la qualité de l’alignement.
rule align_sort_index:
    input:
        bam=get_bam_path,
        ref=REF
    output:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.aligned.sorted.bam.bai",
        flagstat=f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.flagstat.txt",
        alignlog=f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.align.log"
    threads: config["threads"]["align_sort_index"]
    log:
        f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.snakemake.log"
    resources:
        mem_mb=64000,
        runtime=12 * 60,
        tmpdir="/env/cns/proj/projet_DWS/scratch/adjiplif"
    shell:
        r"""
        set -euo pipefail
        module load dorado
        module load samtools

        mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})

        workdir=$(mktemp -d {resources.tmpdir}/${{USER}}.{wildcards.sample}.XXXXXX)

        cp {input.bam} "$workdir/input.bam"
        cp {input.ref} "$workdir/ref.fna"

        dorado aligner "$workdir/ref.fna" "$workdir/input.bam" \
          > "$workdir/{wildcards.sample}.aligned.unsorted.bam" \
          2> "$workdir/{wildcards.sample}.align.log"

        samtools quickcheck -v "$workdir/{wildcards.sample}.aligned.unsorted.bam"

        samtools sort -@ {threads} \
          -o "$workdir/{wildcards.sample}.aligned.sorted.bam" \
          "$workdir/{wildcards.sample}.aligned.unsorted.bam"

        samtools index "$workdir/{wildcards.sample}.aligned.sorted.bam"
        samtools quickcheck -v "$workdir/{wildcards.sample}.aligned.sorted.bam"
        samtools flagstat "$workdir/{wildcards.sample}.aligned.sorted.bam" > "$workdir/{wildcards.sample}.flagstat.txt"

        cp "$workdir/{wildcards.sample}.aligned.sorted.bam" {output.bam}
        cp "$workdir/{wildcards.sample}.aligned.sorted.bam.bai" {output.bai}
        cp "$workdir/{wildcards.sample}.flagstat.txt" {output.flagstat}
        cp "$workdir/{wildcards.sample}.align.log" {output.alignlog}

        rm -rf "$workdir"
        """    


# =============================================================================
# Création des listes de contigs par organisme
# =============================================================================
# À partir de l’index FASTA (.fai), cette règle extrait les identifiants de
# contigs correspondant au puceron nucléaire, à Serratia et à la mitochondrie.
# Ces listes servent ensuite à séparer les BAM par organisme.
rule make_contig_lists_from_new_reference:
    input:
        fai=REF + ".fai"
    output:
        aphid="contig_aphid_list",
        serratia="contig_serratia_list",
        mito="contig_mito_aphid_list"
    log:
        "logs/make_contig_lists_from_new_reference.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p logs

        echo "Création des listes de contigs depuis {input.fai}" > {log}

        # Aphid nucléaire : chromosomes/scaffolds du puceron
        awk '$1 ~ /^NC_042/ || $1 ~ /^NW_/ {{print $1}}' {input.fai} > {output.aphid}

        # Mitochondrie du puceron
        awk '$1 == "NC_011594.1" {{print $1}}' {input.fai} > {output.mito}

        # Serratia : chromosome + plasmides dans ta nouvelle référence
        awk '$1 ~ /^NZ_CP050855/ || $1 ~ /^NZ_CP050856/ || $1 ~ /^NZ_CP050857/ {{print $1}}' {input.fai} > {output.serratia}

        echo "Nombre contigs aphid:" >> {log}
        wc -l {output.aphid} >> {log}

        echo "Nombre contigs mito:" >> {log}
        wc -l {output.mito} >> {log}

        echo "Nombre contigs serratia:" >> {log}
        wc -l {output.serratia} >> {log}

        echo "Premiers contigs aphid:" >> {log}
        head {output.aphid} >> {log}

        echo "Contigs Serratia:" >> {log}
        cat {output.serratia} >> {log}

        echo "Contigs mito:" >> {log}
        cat {output.mito} >> {log}
        """


# =============================================================================
# Séparation des BAM par organisme
# =============================================================================
# Après alignement sur la référence combinée, cette règle extrait les reads
# alignés sur les contigs d’un organisme donné. Elle produit donc un BAM
# spécifique au puceron ou à Serratia.
rule split_bam_by_organism:
    input:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.aligned.sorted.bam.bai",
        contigs=lambda wildcards: CONTIG_LISTS[wildcards.organism]
    output:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam.bai"
    log:
        "results/{sample}/{sample}.{organism}.split_bam.log"
    threads: config["threads"]["split_bam_by_organism"]
    resources:
        mem_mb=16000,
        runtime=4 * 60,
        tmpdir="/env/cns/proj/projet_DWS/scratch/adjiplif"
    shell:
        r"""
        set -euo pipefail
        module load samtools

        mkdir -p $(dirname {output.bam}) results/{wildcards.sample}

        workdir=$(mktemp -d {resources.tmpdir}/${{USER}}.{wildcards.sample}.{wildcards.organism}.split.XXXXXX)

        echo "Contig list: {input.contigs}" > {log}
        echo "First contigs:" >> {log}
        head {input.contigs} >> {log}

        samtools view -b -@ {threads} \
          -o "$workdir/output.bam" \
          {input.bam} \
          $(tr -d '\r' < {input.contigs}) \
          2>> {log}

        samtools index "$workdir/output.bam"
        samtools quickcheck -v "$workdir/output.bam"

        echo "Flagstat:" >> {log}
        samtools flagstat "$workdir/output.bam" >> {log}

        cp "$workdir/output.bam" {output.bam}
        cp "$workdir/output.bam.bai" {output.bai}

        rm -rf "$workdir"
        """
		
# =============================================================================
# Statistiques d’alignement par contig et par organisme
# =============================================================================
# Cette règle utilise samtools idxstats pour obtenir le nombre de reads alignés
# sur chaque contig de l’organisme analysé.
rule idxstats_by_organism:
    input:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam.bai"
    output:
        txt="results/{sample}/{sample}.{organism}.idxstats.txt",
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        r"""
        module load samtools
        samtools idxstats {input.bam} > {output.txt}
        """		
		

# =============================================================================
# Résumé modkit par organisme
# =============================================================================
# Cette règle produit un résumé global des bases modifiées détectées dans le BAM
# d’un organisme donné.
rule modkit_summary_by_organism:
    input:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam.bai"
    output:
        txt="results/{sample}/{sample}.{organism}.modkit_summary.txt"
    threads: config["threads"]["modkit_summary_by_organism"]
    resources:
        mem_mb=16000,
        runtime=4 * 60,
        tmpdir="/env/cns/proj/projet_DWS/scratch/adjiplif"
    shell:
        r"""
        module load modkit
        set -euo pipefail
        workdir=$(mktemp -d {resources.tmpdir}/${{USER}}.{wildcards.sample}.{wildcards.organism}.summary.XXXXXX)

        cp {input.bam} "$workdir/input.bam"
        cp {input.bai} "$workdir/input.bam.bai"

        modkit summary "$workdir/input.bam" > "$workdir/summary.txt"
        cp "$workdir/summary.txt" {output.txt}

        rm -rf "$workdir"
        """




# QC before pileup  on utilise flagstat
# =============================================================================
# QC avant pileup de méthylation
# =============================================================================
# Cette règle résume le taux d’alignement et le nombre de reads alignés à partir
# des fichiers flagstat.
rule plot_qc_before_pileup:
    input:
        expand(f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.flagstat.txt", sample=SAMPLES)
    output:
        "results/qc/qc_before_pileup.png"
    run:
        import re
        import pandas as pd
        import matplotlib.pyplot as plt

        data = []

        for f in input:
            sample = f.split("/")[-2]
            txt = open(f).read()

            m = re.search(r'(\d+) \+ \d+ mapped \(([\d\.]+)%', txt)

            if m:
                mapped = int(m.group(1))
                rate = float(m.group(2))
            else:
                mapped = 0
                rate = 0

            data.append([sample, mapped, rate])

        df = pd.DataFrame(data, columns=["sample", "mapped", "rate"])

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].bar(df["sample"], df["rate"])
        axes[0].set_title("Taux d'alignement (%)")

        axes[1].bar(df["sample"], df["mapped"])
        axes[1].set_title("Reads alignés")

        for ax in axes:
            ax.tick_params(axis="x", rotation=90)

        plt.tight_layout()
        plt.savefig(output[0])

#Couverture par organism
# =============================================================================
# Couverture du génome par organisme
# =============================================================================
# Cette règle calcule la couverture avec samtools coverage, puis conserve
# uniquement les lignes correspondant aux contigs de l’organisme étudié.
rule coverage_by_organism:
    input:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam.bai",
        contig_organism=lambda wildcards: CONTIG_LISTS[wildcards.organism]
    output:
        txt="results/{sample}/{sample}.{organism}.coverage.tsv"
    log:
        "results/{sample}/{sample}.{organism}.coverage.log"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=4 * 60     		
    shell:
        r"""
        set -euo pipefail
        module load samtools

        mkdir -p results/{wildcards.sample}

        samtools coverage {input.bam} > tmp.{wildcards.sample}.{wildcards.organism}.coverage

        grep -Ff {input.contig_organism} tmp.{wildcards.sample}.{wildcards.organism}.coverage > {output.txt}

        rm -f tmp.{wildcards.sample}.{wildcards.organism}.coverage
        """
# =============================================================================
# Figure QC de couverture par organisme
# =============================================================================
# Cette règle agrège les tables de couverture et produit une figure montrant :
#   - le pourcentage du génome couvert ;
#   - la profondeur moyenne de séquençage.
rule plot_qc_coverage_by_organism:
    input:
        expand("results/{sample}/{sample}.{organism}.coverage.tsv",
               sample=SAMPLES, organism=ORGANISMS)
    output:
        "results/qc/qc_coverage_by_organism.png"
    run:
        import pandas as pd, matplotlib.pyplot as plt

        rows = []

        for f in input:
            sample = f.split("/")[-2]
            filename = f.split("/")[-1]
            organism = filename.replace(f"{sample}.", "").replace(".coverage.tsv", "")
            df = pd.read_csv(f, sep=r"\s+")
            # agrégation sur tous les contigs de l’organisme
            total_len = (df["endpos"] - df["startpos"] + 1).sum()
            total_covbases = df["covbases"].sum()
            weighted_meandepth = ((df["endpos"] - df["startpos"] + 1) * df["meandepth"]).sum() / total_len
            total_numreads = df["numreads"].sum()

            pct_cov = 100 * total_covbases / total_len if total_len > 0 else 0

            rows.append([sample, organism, total_numreads, pct_cov, weighted_meandepth])

        out = pd.DataFrame(rows, columns=["sample", "organism", "numreads", "pct_cov", "meandepth"])
        labels = out["sample"] + "." + out["organism"]

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        axes[0].bar(labels, out["pct_cov"])
        axes[0].set_title("% génome couvert par organisme")

        axes[1].bar(labels, out["meandepth"])
        axes[1].set_title("Profondeur moyenne par organisme")

        for ax in axes:
            ax.tick_params(axis="x", rotation=90)

        plt.tight_layout()
        plt.savefig(output[0])


# Retourne le code modkit correspondant à la marque de méthylation demandée.
def mod_base(wildcards):
    return MOD_BASES[wildcards.modtype]
		

# =============================================================================
# Appel des marques de méthylation avec modkit pileup
# =============================================================================
# Cette règle utilise les probabilités de bases modifiées contenues dans les BAM
# Nanopore pour produire un fichier bedMethyl brut.
# Les marques testées sont définies dans MOD_BASES : 4mC, 5mC et 6mA.
rule modkit_pileup:
    input:
        bam=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam",
        bai=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.aligned.sorted.bam.bai",
        ref=REF,
        fai=REF + ".fai"
    output:
        bed=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.raw.bed"
    log:
        "results/{sample}/{sample}.{organism}.{modtype}.modkit_pileup.log"
    threads: config["threads"]["modkit_pileup"]
    resources:
        mem_mb=256000,
        runtime=24 * 60,
        tmpdir="/env/cns/proj/projet_DWS/scratch/adjiplif"
    params:
        extra=MODKIT_EXTRA,
        modbase=mod_base
    shell:
        r"""
        set -euo pipefail
        module load modkit/0.6.2

        workdir=$(mktemp -d {resources.tmpdir}/${{USER}}.{wildcards.sample}.{wildcards.organism}.{wildcards.modtype}.pileup.XXXXXX)

        cp {input.bam} "$workdir/input.bam"
        cp {input.bai} "$workdir/input.bam.bai"
        cp {input.ref} "$workdir/ref.fna"
        cp {input.fai} "$workdir/ref.fna.fai"

        modkit pileup "$workdir/input.bam" "$workdir/output.bed" \
          --ref "$workdir/ref.fna" \
          --modified-bases {params.modbase} \
          {params.extra} \
          > {log} 2>&1

        cp "$workdir/output.bed" {output.bed}

        rm -rf "$workdir"
        """
	


# =============================================================================
# Filtrage des fichiers bedMethyl
# =============================================================================
# Cette règle conserve uniquement les positions ayant une couverture suffisante
# et des comptes cohérents :
#   coverage = reads modifiés + reads canoniques
#   count_other = 0
# Cela évite de garder des positions ambiguës ou de mauvaise qualité.
rule filter_mods_bed_by_coverage:
    input:
        bed=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.raw.bed"
    output:
        bed=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed"
    log:
        "results/{sample}/{sample}.{organism}.{modtype}.filter_cov.log"
    params:
        min_cov=MODKIT_MIN_COV
    shell:
        r"""
        set -euo pipefail

        mkdir -p results/{wildcards.sample}

        echo "Filtrage par couverture minimale + validité bedMethyl" > {log}
        echo "Input: {input.bed}" >> {log}
        echo "Output: {output.bed}" >> {log}
        echo "min_cov: {params.min_cov}" >> {log}
        echo "Condition gardée: coverage >= min_cov ET coverage == count_modified + count_canonical ET count_other == 0" >> {log}

        echo "Nombre de lignes avant filtrage:" >> {log}
        grep -vc "^#" {input.bed} >> {log} || true

        awk -v mincov={params.min_cov} '
          BEGIN {{OFS="\t"}}

          /^#/ {{
            print;
            next;
          }}

          NF >= 14 {{
            cov = $10 + 0;
            n_mod = $12 + 0;
            n_canon = $13 + 0;
            n_other = $14 + 0;

            if (cov >= mincov && cov == (n_mod + n_canon) && n_other == 0) {{
              print;
            }}
          }}
        ' {input.bed} > {output.bed}

        echo "Nombre de lignes après filtrage:" >> {log}
        grep -vc "^#" {output.bed} >> {log} || true

        echo "Nombre de lignes invalides retirées car coverage != modified + canonical ou count_other != 0:" >> {log}
        awk '
          /^#/ {{next}}
          NF >= 14 {{
            cov = $10 + 0;
            n_mod = $12 + 0;
            n_canon = $13 + 0;
            n_other = $14 + 0;
            if (cov != (n_mod + n_canon) || n_other != 0) bad++;
          }}
          END {{print bad + 0}}
        ' {input.bed} >> {log}
        """
 
# =============================================================================
# Compression et indexation des fichiers bedMethyl
# =============================================================================
# Cette règle compresse les fichiers BED avec bgzip et crée un index tabix.
# Ces fichiers indexés sont nécessaires pour les comparaisons différentielles
# avec modkit dmr pair.
rule bgzip_tabix_mods_bed:
    input:
        bed=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed"
    output:
        bed_gz=f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed.gz",
        tbi=f"{SCRATCH_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed.gz.tbi"
    log:
        "results/{sample}/{sample}.{organism}.{modtype}.bgzip_tabix.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=12*60
    shell:
        r"""
        set -euo pipefail
        bgzip -c {input.bed} > {output.bed_gz} 2> {log}
        tabix -f -p bed {output.bed_gz} >> {log} 2>&1
        """

# =============================================================================
# Conversion modkit vers methylKit
# =============================================================================
# methylKit attend une table simple contenant chromosome, position, brin,
# couverture et fréquence de méthylation. Cette règle convertit les fichiers
# bedMethyl filtrés vers ce format.
# Cas particulier : pour la 6mA chez le puceron, les sites sont agrégés par
# fenêtres afin de réduire le bruit lié au faible niveau de signal.
rule convert_modkit_to_methylkit:
    input:
        bed=f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed"
    output:
        tsv=f"{METHYLKIT_DIR}/input/{{sample}}.{{organism}}.{{modtype}}.methylkit.tsv"
    params:
        mincov=config["modkit"].get("min_coverage", 5),
        win=10000,
        mode=lambda wc: "window" if (wc.organism == "aphid" and wc.modtype == "6mA") else "base"
    log:
        f"{METHYLKIT_DIR}/logs/input/{{sample}}.{{organism}}.{{modtype}}.convert.log"
    shell:
        r"""
        mkdir -p $(dirname {output.tsv}) $(dirname {log})

        echo "Conversion modkit -> methylKit" > {log}
        echo "sample={wildcards.sample}" >> {log}
        echo "organism={wildcards.organism}" >> {log}
        echo "modtype={wildcards.modtype}" >> {log}
        echo "mode={params.mode}" >> {log}
        echo "mincov={params.mincov}" >> {log}
        echo "window={params.win}" >> {log}

        if [ "{params.mode}" = "window" ]; then

            echo "Mode spécial aphid 6mA : agrégation par fenêtres" >> {log}

            awk -v OFS="\t" -v MINCOV="{params.mincov}" -v W="{params.win}" '
            BEGIN {{
                print "chr","start","end","strand","coverage","freqC"
            }}
            NR > 1 {{
                chr=$1
                pos=$2
                cov=$10
                freq=$11
                n_mod=$12
                n_canon=$13
                n_other=$14

                if (cov >= MINCOV && n_other == 0 && cov == (n_mod + n_canon)) {{
                    wstart = int(pos / W) * W + 1
                    wend = wstart + W - 1
                    key = chr OFS wstart OFS wend OFS "+"

                    cov_sum[key] += cov
                    mod_sum[key] += n_mod
                }}
            }}
            END {{
                for (key in cov_sum) {{
                    if (cov_sum[key] > 0) {{
                        freq_win = 100 * mod_sum[key] / cov_sum[key]
                        print key, cov_sum[key], freq_win
                    }}
                }}
            }}
            ' {input.bed} \
            | awk 'NR==1 {{print; next}} {{print | "sort -k1,1 -k2,2n"}}' \
            > {output.tsv}

        else

            echo "Mode normal : base par base" >> {log}

            awk -v OFS="\t" -v MINCOV="{params.mincov}" '
            BEGIN {{
                print "chr","start","end","strand","coverage","freqC"
            }}
            NR > 1 {{
                chr=$1
                start=$2 + 1
                end=$3
                strand=$6
                cov=$10
                freq=$11
                n_mod=$12
                n_canon=$13
                n_other=$14

                if (strand != "+" && strand != "-") {{
                    strand = "+"
                }}

                if (cov >= MINCOV && n_other == 0 && cov == (n_mod + n_canon)) {{
                    print chr,start,end,strand,cov,freq
                }}
            }}
            ' {input.bed} > {output.tsv}

        fi

        echo "Nombre de lignes produites :" >> {log}
        wc -l {output.tsv} >> {log}
        """

# ============================================================
# Figure methylKit : nombre de DMR hyper/hypométhylées
# Version classique avec échelles libres
#
# Entrées :
#   results/methylkit/diff/{comp}.{organism}.{modtype}.counts.tsv
#
# Sorties :
#   results/methylkit/tables/methylkit_hyper_hypo_counts_all.tsv
#   results/methylkit/plots/methylkit_hyper_hypo_all_comparisons.png
# ============================================================

# =============================================================================
# Figure methylKit : DMR hyper- et hypométhylées par comparaison
# =============================================================================
# Cette règle lit les fichiers counts.tsv produits après l’analyse methylKit et
# génère un histogramme avec échelles libres entre les panneaux.
rule plot_methylkit_hyper_hypo_all_comparisons:
    input:
        counts=METHYLKIT_COUNTS
    output:
        table=f"{METHYLKIT_DIR}/tables/methylkit_hyper_hypo_counts_all.tsv",
        plot=f"{METHYLKIT_DIR}/plots/methylkit_hyper_hypo_all_comparisons.png"
    log:
        f"{METHYLKIT_DIR}/logs/methylkit_hyper_hypo_all_comparisons.log"
    resources:
        mem_mb=8000,
        runtime=60
    run:
        import os
        import subprocess

        os.makedirs(f"{METHYLKIT_DIR}/tables", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/plots", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/logs", exist_ok=True)

        r_script = f"{METHYLKIT_DIR}/logs/plot_methylkit_hyper_hypo_all_comparisons.R"
        counts_files = " ".join(list(input.counts))

        r_code = r'''
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

files <- strsplit("__INPUT_COUNTS__", " ")[[1]]

message("Nombre de fichiers methylKit counts reçus : ", length(files))

all_counts <- list()

for (f in files) {

    message("Lecture : ", f)

    if (!file.exists(f) || file.info(f)$size == 0) {
        message("Fichier absent ou vide, ignoré : ", f)
        next
    }

    df <- tryCatch(
        read_tsv(f, show_col_types = FALSE),
        error = function(e) {
            message("Erreur lecture : ", f)
            message(e$message)
            return(NULL)
        }
    )

    if (is.null(df) || nrow(df) == 0) {
        message("Dataframe vide, ignoré : ", f)
        next
    }

    required_cols <- c("comparison", "organism", "modtype", "n_hyper", "n_hypo")
    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
        stop(
            "Colonnes manquantes dans : ", f,
            "\nColonnes manquantes : ", paste(missing_cols, collapse = ", "),
            "\nColonnes disponibles : ", paste(colnames(df), collapse = ", ")
        )
    }

    df2 <- df %>%
        select(comparison, organism, modtype, n_hyper, n_hypo) %>%
        mutate(
            n_hyper = as.numeric(n_hyper),
            n_hypo = as.numeric(n_hypo),
            total = n_hyper + n_hypo
        )

    all_counts[[length(all_counts) + 1]] <- df2
}

res <- bind_rows(all_counts)

if (nrow(res) == 0) {
    stop("Aucune donnée utilisable pour produire le graphique.")
}

write_tsv(res, "__OUTPUT_TABLE__")

plot_df <- res %>%
    pivot_longer(
        cols = c("n_hyper", "n_hypo"),
        names_to = "direction",
        values_to = "n_dmr"
    ) %>%
    mutate(
        direction = recode(
            direction,
            "n_hyper" = "Hyperméthylé",
            "n_hypo" = "Hypométhylé"
        ),
        organism = recode(
            organism,
            "aphid" = "Acyrthosiphon pisum",
            "serratia" = "Serratia symbiotica"
        ),
        comparison = factor(comparison, levels = unique(comparison)),
        modtype = factor(modtype, levels = c("4mC", "5mC", "6mA")),
        direction = factor(direction, levels = c("Hyperméthylé", "Hypométhylé"))
    )

p <- ggplot(plot_df, aes(x = comparison, y = n_dmr, fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(modtype ~ organism, scales = "free_y") +
    labs(
        title = "Nombre de DMR hyper- et hypométhylées par comparaison",
        subtitle = "Détection basée sur methylKit",
        x = "Comparaison",
        y = "Nombre de DMR retenues",
        fill = "Direction"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 7),
        strip.text = element_text(face = "plain"),
        legend.position = "bottom"
    )

ggsave(
    filename = "__OUTPUT_PLOT__",
    plot = p,
    width = 18,
    height = 9,
    dpi = 300
)

message("Table écrite : __OUTPUT_TABLE__")
message("Graphique écrit : __OUTPUT_PLOT__")
'''

        r_code = r_code.replace("__INPUT_COUNTS__", counts_files)
        r_code = r_code.replace("__OUTPUT_TABLE__", str(output.table))
        r_code = r_code.replace("__OUTPUT_PLOT__", str(output.plot))

        with open(r_script, "w") as f:
            f.write(r_code)

        cmd = f"""
        set -euo pipefail
        module load bioconductor 2>/dev/null || true
        Rscript {r_script} > {log[0]} 2>&1
        """

        subprocess.check_call(cmd, shell=True, executable="/bin/bash")
## fonctions d'entrée pour mon modkit_dmr pair
# Liste des fichiers bedMethyl compressés du groupe A pour modkit dmr pair.
def pair_a_beds(wildcards):
    group_a, _ = get_samples_for_comparison(wildcards.comp)
    return [
        f"{SCRATCH_RESULTS}/{s}/{s}.{wildcards.organism}.{wildcards.modtype}.mods.bed.gz"
        for s in group_a
    ]

# Liste des fichiers bedMethyl compressés du groupe B pour modkit dmr pair.
def pair_b_beds(wildcards):
    _, group_b = get_samples_for_comparison(wildcards.comp)
    return [
        f"{SCRATCH_RESULTS}/{s}/{s}.{wildcards.organism}.{wildcards.modtype}.mods.bed.gz"
        for s in group_b
    ]

def pair_a_tbis(wildcards):
    group_a, _ = get_samples_for_comparison(wildcards.comp)
    return [
        f"{SCRATCH_RESULTS}/{s}/{s}.{wildcards.organism}.{wildcards.modtype}.mods.bed.gz.tbi"
        for s in group_a
    ]

def pair_b_tbis(wildcards):
    _, group_b = get_samples_for_comparison(wildcards.comp)
    return [
        f"{SCRATCH_RESULTS}/{s}/{s}.{wildcards.organism}.{wildcards.modtype}.mods.bed.gz.tbi"
        for s in group_b
    ]

def pair_args_a(wildcards):
    group_a, _ = get_samples_for_comparison(wildcards.comp)
    return " ".join(
        f"-a {SCRATCH_RESULTS}/{s}/{s}.{wildcards.organism}.{wildcards.modtype}.mods.bed.gz"
        for s in group_a
    )

def pair_args_b(wildcards):
    _, group_b = get_samples_for_comparison(wildcards.comp)
    return " ".join(
        f"-b {SCRATCH_RESULTS}/{s}/{s}.{wildcards.organism}.{wildcards.modtype}.mods.bed.gz"
        for s in group_b
    )


# Retourne la base canonique analysée par modkit dmr pair :
# C pour les cytosines méthylées, A pour la 6mA.
def dmr_base(wildcards):
    return DMR_BASE[wildcards.modtype]

# Pour la 4mC, précise à modkit le code de modification non standard utilisé.
def dmr_assign_code(wildcards):
    if wildcards.modtype == "4mC":
        return "--assign-code 21839:C"
    return ""



## ##QC avant modkit_dmr on utilise mods.bed pour évaluer la qualité de la méthylation
# =============================================================================
# QC avant détection des DMR
# =============================================================================
# Cette règle vérifie les fichiers bedMethyl filtrés avant les analyses
# différentielles : couverture médiane, méthylation moyenne et nombre de sites.
rule plot_qc_before_dmr:
    input:
        expand(
            f"{BIGTMP_RESULTS}/{{sample}}/{{sample}}.{{organism}}.{{modtype}}.mods.bed",
            sample=SAMPLES,
            organism=ORGANISMS,
            modtype=MOD_TYPES
        )
    output:
        "results/qc/qc_before_dmr.png"
    log:
        "results/qc/qc_before_dmr.log"
    run:
        import os
        import pandas as pd
        import matplotlib.pyplot as plt

        os.makedirs("results/qc", exist_ok=True)

        data = []

        with open(log[0], "w") as lg:
            lg.write("QC before DMR\n")
            lg.write("====================\n\n")

            for f in input:
                parts = f.split("/")
                sample = parts[-2]
                filename = parts[-1]

                core = filename.replace(".mods.bed", "")
                prefix = f"{sample}."

                if not core.startswith(prefix):
                    lg.write(f"SKIP nom inattendu: {f}\n")
                    continue

                rest = core[len(prefix):]
                organism, modtype = rest.rsplit(".", 1)

                lg.write(f"\nFichier: {f}\n")
                lg.write(f"Sample={sample}, organism={organism}, modtype={modtype}\n")

                if not os.path.exists(f):
                    lg.write("SKIP: fichier absent\n")
                    data.append([sample, organism, modtype, 0, 0, 0])
                    continue

                if os.path.getsize(f) == 0:
                    lg.write("SKIP: fichier vide\n")
                    data.append([sample, organism, modtype, 0, 0, 0])
                    continue

                try:
                    df = pd.read_csv(f, sep="\t", comment="#", header=None)
                except Exception as e:
                    lg.write(f"SKIP: erreur lecture pandas: {e}\n")
                    data.append([sample, organism, modtype, 0, 0, 0])
                    continue

                if df.empty:
                    lg.write("SKIP: dataframe vide\n")
                    data.append([sample, organism, modtype, 0, 0, 0])
                    continue

                if df.shape[1] < 11:
                    lg.write(f"SKIP: pas assez de colonnes: {df.shape[1]}\n")
                    data.append([sample, organism, modtype, 0, 0, len(df)])
                    continue

                try:
                    cov = pd.to_numeric(df.iloc[:, 9], errors="coerce").dropna()
                    meth = pd.to_numeric(df.iloc[:, 10], errors="coerce").dropna()

                    median_cov = cov.median() if len(cov) > 0 else 0
                    mean_meth = meth.mean() if len(meth) > 0 else 0
                    n_sites = len(df)

                    lg.write(f"OK: n_sites={n_sites}, median_cov={median_cov}, mean_meth={mean_meth}\n")

                    data.append([
                        sample,
                        organism,
                        modtype,
                        median_cov,
                        mean_meth,
                        n_sites
                    ])

                except Exception as e:
                    lg.write(f"SKIP: erreur calcul: {e}\n")
                    data.append([sample, organism, modtype, 0, 0, len(df)])

        out = pd.DataFrame(
            data,
            columns=[
                "sample",
                "organism",
                "modtype",
                "median_cov",
                "mean_meth",
                "n_sites"
            ]
        )

        out.to_csv("results/qc/qc_before_dmr_table.tsv", sep="\t", index=False)

        if out.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "Aucune donnée disponible", ha="center", va="center")
            ax.axis("off")
            plt.savefig(output[0], dpi=300)
        else:
            labels = out["sample"] + "." + out["organism"] + "." + out["modtype"]

            fig, axes = plt.subplots(1, 3, figsize=(20, 6))

            axes[0].bar(labels, out["median_cov"])
            axes[0].set_title("Couverture médiane")
            axes[0].set_ylabel("Couverture")

            axes[1].bar(labels, out["mean_meth"])
            axes[1].set_title("Méthylation moyenne")
            axes[1].set_ylabel("% méthylation")

            axes[2].bar(labels, out["n_sites"])
            axes[2].set_title("Nombre de sites")
            axes[2].set_ylabel("n sites")

            for ax in axes:
                ax.tick_params(axis="x", rotation=90)

            plt.tight_layout()
            plt.savefig(output[0], dpi=300)


# =============================================================================
# Référence FASTA spécifique à chaque organisme
# =============================================================================
# modkit dmr pair utilise une référence correspondant à l’organisme testé.
# Cette règle extrait donc les contigs de l’organisme depuis la référence combinée.
rule make_organism_ref:
    input:
        ref=REF,
        contigs=lambda wildcards: CONTIG_LISTS[wildcards.organism]
    output:
        fa=f"{SCRATCH}/refs/{{organism}}.only.fa",
        fai=f"{SCRATCH}/refs/{{organism}}.only.fa.fai"
    shell:
        r"""
        set -euo pipefail
        module load samtools
        mkdir -p $(dirname {output.fa})

        samtools faidx {input.ref} $(tr -d '\r' < {input.contigs}) > {output.fa}
        samtools faidx {output.fa}
        """
		
# =============================================================================
# Détection de DMR avec modkit dmr pair
# =============================================================================
# Cette règle compare les groupes A et B définis dans le fichier comparaisons.csv.
# Elle utilise les fichiers bedMethyl compressés/indexés de chaque groupe et
# produit un fichier BED de méthylation différentielle.
rule modkit_dmr_pair:
    input:
        a_beds=pair_a_beds,
        b_beds=pair_b_beds,
        a_tbis=pair_a_tbis,
        b_tbis=pair_b_tbis,
        ref=organism_ref
    output:
        dmr="results/dmr_pair/{comp}.{organism}.{modtype}.bed"
    log:
        shell="results/dmr_pair/{comp}.{organism}.{modtype}.shell.log",
        modkit="results/dmr_pair/{comp}.{organism}.{modtype}.modkit.log"
    threads: config["threads"]["modkit_dmr_pair"]
    resources:
        mem_mb=128000,
        runtime=24 * 60
    params:
        args_a=pair_args_a,
        args_b=pair_args_b,
        base=dmr_base,
        assign_code=dmr_assign_code,
        tmpdir="/env/cns/proj/projet_DWS/scratch/adjiplif"
    shell:
        r"""
        set -euo pipefail
        module load modkit/0.6.2

        mkdir -p results/dmr_pair
        mkdir -p {params.tmpdir}

        workdir=$(mktemp -d {params.tmpdir}/${{USER}}.{wildcards.comp}.{wildcards.organism}.{wildcards.modtype}.dmr.XXXXXX)

        export TMPDIR="$workdir"
        export TMP="$workdir"
        export TEMP="$workdir"

        echo "DMR pair: {wildcards.comp} {wildcards.organism} {wildcards.modtype}" > {log.shell}
        echo "TMPDIR=$TMPDIR" >> {log.shell}
        echo "Output: {output.dmr}" >> {log.shell}
        echo "Reference: {input.ref}" >> {log.shell}
        echo "Base: {params.base}" >> {log.shell}
        echo "Assign code: {params.assign_code}" >> {log.shell}
        echo "A inputs: {input.a_beds}" >> {log.shell}
        echo "B inputs: {input.b_beds}" >> {log.shell}

        echo "MODKIT VERSION:" >> {log.shell}
        which modkit >> {log.shell} 2>&1
        modkit --version >> {log.shell} 2>&1

        modkit dmr pair \
          {params.args_a} \
          {params.args_b} \
          -o {output.dmr} \
          --ref {input.ref} \
          --base {params.base} \
          {params.assign_code} \
          --threads {threads} \
          --log-filepath {log.modkit} \
          >> {log.shell} 2>&1

        rm -rf "$workdir"
        """
# =============================================================================
# Figure methylKit avec échelle fixée sur la 6mA
# =============================================================================
# Cette version reprend les counts methylKit mais fixe l’axe Y selon l’amplitude
# observée en 6mA. Elle permet de mieux visualiser les faibles nombres de DMR.
rule plot_methylkit_hyper_hypo_scale_6mA:
    input:
        counts=METHYLKIT_COUNTS
    output:
        table=f"{METHYLKIT_DIR}/tables/methylkit_hyper_hypo_counts_all_scale_6mA.tsv",
        plot=f"{METHYLKIT_DIR}/plots/methylkit_hyper_hypo_all_comparisons_scale_6mA.png"
    log:
        f"{METHYLKIT_DIR}/logs/methylkit_hyper_hypo_scale_6mA.log"
    run:
        import os
        import subprocess

        os.makedirs(f"{METHYLKIT_DIR}/tables", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/plots", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/logs", exist_ok=True)

        r_script = f"{METHYLKIT_DIR}/logs/plot_methylkit_hyper_hypo_scale_6mA.R"
        counts_files = " ".join(list(input.counts))

        r_code = r'''
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

files <- strsplit("__INPUT_COUNTS__", " ")[[1]]

message("Nombre de fichiers methylKit counts reçus : ", length(files))

all_counts <- list()

for (f in files) {

    message("Lecture : ", f)

    if (!file.exists(f) || file.info(f)$size == 0) {
        message("Fichier absent ou vide, ignoré : ", f)
        next
    }

    df <- tryCatch(
        read_tsv(f, show_col_types = FALSE),
        error = function(e) {
            message("Erreur lecture : ", f)
            return(NULL)
        }
    )

    if (is.null(df) || nrow(df) == 0) {
        message("Dataframe vide, ignoré : ", f)
        next
    }

    required_cols <- c("comparison", "organism", "modtype", "n_hyper", "n_hypo")
    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
        stop(
            "Colonnes manquantes dans : ", f,
            "\nColonnes manquantes : ", paste(missing_cols, collapse = ", "),
            "\nColonnes disponibles : ", paste(colnames(df), collapse = ", ")
        )
    }

    df2 <- df %>%
        select(comparison, organism, modtype, n_hyper, n_hypo) %>%
        mutate(
            n_hyper = as.numeric(n_hyper),
            n_hypo = as.numeric(n_hypo),
            total = n_hyper + n_hypo
        )

    all_counts[[length(all_counts) + 1]] <- df2
}

res <- bind_rows(all_counts)

if (nrow(res) == 0) {
    stop("Aucune donnée utilisable pour produire le graphique.")
}

write_tsv(res, "__OUTPUT_TABLE__")

plot_df <- res %>%
    pivot_longer(
        cols = c("n_hyper", "n_hypo"),
        names_to = "direction",
        values_to = "n_dmr"
    ) %>%
    mutate(
        direction = recode(
            direction,
            "n_hyper" = "Hyperméthylé",
            "n_hypo" = "Hypométhylé"
        ),
        comparison = factor(comparison, levels = unique(comparison)),
        modtype = factor(modtype, levels = c("4mC", "5mC", "6mA")),
        direction = factor(direction, levels = c("Hyperméthylé", "Hypométhylé"))
    )

max_6mA <- plot_df %>%
    filter(modtype == "6mA") %>%
    summarise(max_val = max(n_dmr, na.rm = TRUE)) %>%
    pull(max_val)

if (length(max_6mA) == 0 || is.na(max_6mA) || max_6mA <= 0) {
    message("Aucun DMR 6mA détecté. Échelle fixée par défaut à 1.")
    max_6mA <- 1
}

message("Échelle Y fixée sur 6mA : 0 à ", max_6mA)

p <- ggplot(plot_df, aes(x = comparison, y = n_dmr, fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(modtype ~ organism, scales = "fixed") +
    coord_cartesian(ylim = c(0, max_6mA)) +
    labs(
        title = "Nombre de DMR hyper- et hypométhylées par comparaison",
        subtitle = "Détection basée sur methylKit — échelle fixée sur 6mA",
        x = "Comparaison",
        y = "Nombre de DMR retenues",
        fill = "Direction"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 7),
        strip.text = element_text(face = "plain"),
        legend.position = "bottom"
    )

ggsave(
    filename = "__OUTPUT_PLOT__",
    plot = p,
    width = 18,
    height = 9,
    dpi = 300
)

message("Table écrite : __OUTPUT_TABLE__")
message("Graphique écrit : __OUTPUT_PLOT__")
'''

        r_code = r_code.replace("__INPUT_COUNTS__", counts_files)
        r_code = r_code.replace("__OUTPUT_TABLE__", str(output.table))
        r_code = r_code.replace("__OUTPUT_PLOT__", str(output.plot))

        with open(r_script, "w") as f:
            f.write(r_code)

        cmd = f"""
        set -euo pipefail
        module load bioconductor 2>/dev/null || true
        Rscript {r_script} > {log[0]} 2>&1
        """

        subprocess.check_call(cmd, shell=True, executable="/bin/bash")
		
# =============================================================================
# Annotation des DMR methylKit aux gènes
# =============================================================================
# Cette règle lit les DMR hyper/hypo produites par methylKit, puis les associe
# aux gènes voisins ou chevauchants à l’aide d’un fichier BED de gènes.
# Elle produit une table annotée et une liste unique de gènes associés aux DMR.
rule annotate_methylkit_dmr:
    input:
        hyper=f"{METHYLKIT_DIR}/diff/{{comp}}.{{organism}}.{{modtype}}.hyper.tsv",
        hypo=f"{METHYLKIT_DIR}/diff/{{comp}}.{{organism}}.{{modtype}}.hypo.tsv",
        genes=lambda wc: GENE_BEDS[wc.organism]
    output:
        annotated=f"{METHYLKIT_DIR}/annotation/{{comp}}.{{organism}}.{{modtype}}.dmr_annotated.tsv",
        genes=f"{METHYLKIT_DIR}/annotation/{{comp}}.{{organism}}.{{modtype}}.dmr_genes.txt"
    log:
        f"{METHYLKIT_DIR}/logs/annotation/{{comp}}.{{organism}}.{{modtype}}.annotate_dmr.log"
    params:
        window=DMR_GENE_WINDOW
    run:
        import os
        import subprocess

        os.makedirs(f"{METHYLKIT_DIR}/annotation", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/logs/annotation", exist_ok=True)

        r_script = f"{METHYLKIT_DIR}/logs/annotation/{wildcards.comp}.{wildcards.organism}.{wildcards.modtype}.annotate_dmr.R"

        r_code = r'''
library(readr)
library(dplyr)
library(stringr)

hyper_file <- "__HYPER_FILE__"
hypo_file <- "__HYPO_FILE__"
gene_bed <- "__GENE_BED__"
out_annot <- "__OUT_ANNOT__"
out_genes <- "__OUT_GENES__"
window <- as.numeric("__WINDOW__")

message("Fichier hyper : ", hyper_file)
message("Fichier hypo  : ", hypo_file)
message("Fichier gènes : ", gene_bed)
message("Fenêtre autour des gènes : ", window)

read_dmr_file <- function(file, direction) {

    if (!file.exists(file)) {
        message("Fichier absent : ", file)
        return(tibble())
    }

    if (file.info(file)$size == 0) {
        message("Fichier vide : ", file)
        return(tibble())
    }

    df <- tryCatch(
        read_tsv(file, show_col_types = FALSE),
        error = function(e) {
            message("Erreur lecture : ", file)
            return(tibble())
        }
    )

    if (nrow(df) == 0) {
        message("Aucune ligne dans : ", file)
        return(tibble())
    }

    df$direction <- direction
    return(df)
}

hyper <- read_dmr_file(hyper_file, "hyper")
hypo <- read_dmr_file(hypo_file, "hypo")

dmr <- bind_rows(hyper, hypo)

if (nrow(dmr) == 0) {
    message("Aucune DMR hyper ou hypo à annoter.")
    write_tsv(tibble(), out_annot)
    writeLines(character(0), out_genes)
    quit(save = "no", status = 0)
}

message("Nombre de DMR à annoter : ", nrow(dmr))

if (!all(c("chr", "start", "end") %in% colnames(dmr))) {
    stop(
        "Colonnes attendues absentes dans les fichiers DMR. Il faut chr/start/end.\n",
        "Colonnes disponibles : ", paste(colnames(dmr), collapse = ", ")
    )
}

dmr2 <- dmr %>%
    mutate(
        dmr_chr = as.character(chr),
        dmr_start = as.numeric(start),
        dmr_end = as.numeric(end)
    ) %>%
    filter(!is.na(dmr_chr), !is.na(dmr_start), !is.na(dmr_end))

genes <- read_tsv(
    gene_bed,
    col_names = FALSE,
    show_col_types = FALSE
)

if (ncol(genes) < 4) {
    stop("Le fichier BED des gènes doit avoir au moins 4 colonnes : chr, start, end, attributs")
}

extract_gene_name <- function(x) {

    name <- str_match(x, "Name=([^;]+)")[,2]
    locus <- str_match(x, "locus_tag=([^;]+)")[,2]
    gene <- str_match(x, "gene=([^;]+)")[,2]
    id <- str_match(x, "ID=gene-([^;]+)")[,2]

    out <- ifelse(!is.na(name), name,
           ifelse(!is.na(locus), locus,
           ifelse(!is.na(gene), gene,
           ifelse(!is.na(id), id, x))))

    return(out)
}

genes <- genes %>%
    transmute(
        gene_chr = as.character(X1),
        gene_start = as.numeric(X2),
        gene_end = as.numeric(X3),
        gene_info = as.character(X4),
        gene_id = extract_gene_name(as.character(X4)),
        gene_strand = ifelse(ncol(genes) >= 6, as.character(X6), NA_character_)
    ) %>%
    filter(!is.na(gene_chr), !is.na(gene_start), !is.na(gene_end), !is.na(gene_id))

message("Nombre de gènes dans le BED : ", nrow(genes))

annot_list <- list()

for (i in seq_len(nrow(dmr2))) {

    one <- dmr2[i, ]

    hits <- genes %>%
        filter(gene_chr == one$dmr_chr) %>%
        filter(gene_end >= (one$dmr_start - window)) %>%
        filter(gene_start <= (one$dmr_end + window))

    if (nrow(hits) == 0) {

        tmp <- one %>%
            mutate(
                gene_id = NA_character_,
                gene_info = NA_character_,
                gene_chr = NA_character_,
                gene_start = NA_real_,
                gene_end = NA_real_,
                gene_strand = NA_character_,
                distance_to_gene = NA_real_,
                annotation = "intergenic"
            )

    } else {

        tmp <- hits %>%
            mutate(
                distance_to_gene = case_when(
                    gene_end < one$dmr_start ~ one$dmr_start - gene_end,
                    gene_start > one$dmr_end ~ gene_start - one$dmr_end,
                    TRUE ~ 0
                ),
                annotation = ifelse(distance_to_gene == 0, "overlap_gene", "near_gene")
            ) %>%
            bind_cols(one[rep(1, nrow(.)), ])
    }

    annot_list[[i]] <- tmp
}

annot <- bind_rows(annot_list)

write_tsv(annot, out_annot)

gene_ids <- annot %>%
    filter(!is.na(gene_id)) %>%
    pull(gene_id) %>%
    unique() %>%
    sort()

writeLines(gene_ids, out_genes)

message("Nombre de lignes annotées : ", nrow(annot))
message("Nombre de gènes associés : ", length(gene_ids))
message("Table annotée écrite : ", out_annot)
message("Liste de gènes écrite : ", out_genes)
'''

        r_code = r_code.replace("__HYPER_FILE__", str(input.hyper))
        r_code = r_code.replace("__HYPO_FILE__", str(input.hypo))
        r_code = r_code.replace("__GENE_BED__", str(input.genes))
        r_code = r_code.replace("__OUT_ANNOT__", str(output.annotated))
        r_code = r_code.replace("__OUT_GENES__", str(output.genes))
        r_code = r_code.replace("__WINDOW__", str(params.window))

        with open(r_script, "w") as f:
            f.write(r_code)

        cmd = f"""
        set -euo pipefail
        module load bioconductor 2>/dev/null || true
        Rscript {r_script} > {log[0]} 2>&1
        """

        subprocess.check_call(cmd, shell=True, executable="/bin/bash")
		
# =============================================================================
# Enrichissement GO des gènes associés aux DMR methylKit
# =============================================================================
# Cette règle utilise la liste des gènes associés aux DMR et une table gene2GO
# pour rechercher les termes GO surreprésentés.
rule methylkit_go_enrichment:
    input:
        genes=f"{METHYLKIT_DIR}/annotation/{{comp}}.aphid.{{modtype}}.dmr_genes.txt",
        gene2go="results/annotation/aphid_gene2go.tsv"
    output:
        enrichment=f"{METHYLKIT_DIR}/enrichment/{{comp}}.aphid.{{modtype}}.GO_enrichment.tsv"
    log:
        f"{METHYLKIT_DIR}/logs/enrichment/{{comp}}.aphid.{{modtype}}.GO_enrichment.log"
    run:
        import os
        import subprocess

        os.makedirs(f"{METHYLKIT_DIR}/enrichment", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/logs/enrichment", exist_ok=True)

        r_script = f"{METHYLKIT_DIR}/logs/enrichment/{wildcards.comp}.aphid.{wildcards.modtype}.GO_enrichment.R"

        r_code = r'''
library(readr)
library(dplyr)

genes_file <- "__GENES_FILE__"
gene2go_file <- "__GENE2GO_FILE__"
out_file <- "__OUT_FILE__"

message("Liste de gènes DMR : ", genes_file)
message("Fichier gene2go : ", gene2go_file)

if (!file.exists(genes_file)) {
    stop("Fichier de gènes DMR absent : ", genes_file)
}

if (!file.exists(gene2go_file)) {
    stop("Fichier gene2go absent : ", gene2go_file)
}

dmr_genes <- read_lines(genes_file)
dmr_genes <- unique(dmr_genes[dmr_genes != ""])

gene2go <- read_tsv(gene2go_file, show_col_types = FALSE)

if (!all(c("gene_id", "GO_ID") %in% colnames(gene2go))) {
    stop("Le fichier gene2go doit contenir les colonnes gene_id et GO_ID")
}

gene2go <- gene2go %>%
    filter(!is.na(gene_id), !is.na(GO_ID)) %>%
    distinct(gene_id, GO_ID)

universe_genes <- unique(gene2go$gene_id)

dmr_genes_with_go <- intersect(dmr_genes, universe_genes)

message("Nombre de gènes DMR total : ", length(dmr_genes))
message("Nombre de gènes DMR avec GO : ", length(dmr_genes_with_go))
message("Nombre de gènes dans l'univers GO : ", length(universe_genes))

if (length(dmr_genes_with_go) == 0) {
    message("Aucun gène DMR avec annotation GO. Fichier vide produit.")
    write_tsv(
        tibble(
            GO_ID = character(),
            n_dmr_genes_with_GO = integer(),
            n_universe_genes_with_GO = integer(),
            n_dmr_genes_tested = integer(),
            n_universe_genes = integer(),
            pvalue = numeric(),
            padj = numeric(),
            genes = character()
        ),
        out_file
    )
    quit(save = "no", status = 0)
}

go_terms <- unique(gene2go$GO_ID)

res_list <- list()

for (go in go_terms) {

    genes_with_go <- unique(gene2go$gene_id[gene2go$GO_ID == go])

    a <- length(intersect(dmr_genes_with_go, genes_with_go))
    b <- length(setdiff(dmr_genes_with_go, genes_with_go))
    c <- length(setdiff(genes_with_go, dmr_genes_with_go))
    d <- length(setdiff(universe_genes, union(dmr_genes_with_go, genes_with_go)))

    if (a == 0) {
        next
    }

    mat <- matrix(c(a, b, c, d), nrow = 2)

    ft <- fisher.test(mat, alternative = "greater")

    res_list[[length(res_list) + 1]] <- tibble(
        GO_ID = go,
        n_dmr_genes_with_GO = a,
        n_universe_genes_with_GO = length(genes_with_go),
        n_dmr_genes_tested = length(dmr_genes_with_go),
        n_universe_genes = length(universe_genes),
        pvalue = ft$p.value,
        genes = paste(intersect(dmr_genes_with_go, genes_with_go), collapse = ";")
    )
}

res <- bind_rows(res_list)

if (nrow(res) == 0) {
    message("Aucun terme GO avec au moins un gène DMR.")
    write_tsv(
        tibble(
            GO_ID = character(),
            n_dmr_genes_with_GO = integer(),
            n_universe_genes_with_GO = integer(),
            n_dmr_genes_tested = integer(),
            n_universe_genes = integer(),
            pvalue = numeric(),
            padj = numeric(),
            genes = character()
        ),
        out_file
    )
    quit(save = "no", status = 0)
}

res <- res %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    arrange(padj, pvalue)

write_tsv(res, out_file)

message("Nombre de termes GO testés avec au moins un gène DMR : ", nrow(res))
message("Résultat écrit : ", out_file)
'''

        r_code = r_code.replace("__GENES_FILE__", str(input.genes))
        r_code = r_code.replace("__GENE2GO_FILE__", str(input.gene2go))
        r_code = r_code.replace("__OUT_FILE__", str(output.enrichment))

        with open(r_script, "w") as f:
            f.write(r_code)

        cmd = f"""
        set -euo pipefail
        module load bioconductor 2>/dev/null || true
        Rscript {r_script} > {log[0]} 2>&1
        """

        subprocess.check_call(cmd, shell=True, executable="/bin/bash")


# =============================================================================
# Figure d’enrichissement GO
# =============================================================================
# Cette règle transforme les résultats d’enrichissement GO en figure lisible,
# généralement sous forme de barplot des termes les plus significatifs.
rule plot_go_enrichment_barplot:
    input:
        enrichment=f"{METHYLKIT_DIR}/enrichment/{{comp}}.aphid.{{modtype}}.GO_enrichment.tsv"
    output:
        plot=f"{METHYLKIT_DIR}/enrichment/{{comp}}.aphid.{{modtype}}.GO_enrichment_barplot.png"
    log:
        f"{METHYLKIT_DIR}/logs/enrichment/{{comp}}.aphid.{{modtype}}.GO_enrichment_barplot.log"
    params:
        top_n=30
    run:
        import os
        import subprocess

        os.makedirs(f"{METHYLKIT_DIR}/enrichment", exist_ok=True)
        os.makedirs(f"{METHYLKIT_DIR}/logs/enrichment", exist_ok=True)

        r_script = f"{METHYLKIT_DIR}/logs/enrichment/{wildcards.comp}.aphid.{wildcards.modtype}.GO_enrichment_barplot.R"

        r_code = r'''
library(readr)
library(dplyr)
library(ggplot2)

infile <- "__INPUT__"
outfile <- "__OUTPUT__"
top_n <- as.numeric("__TOP_N__")

message("Fichier enrichissement : ", infile)
message("Figure : ", outfile)

df <- read_tsv(infile, show_col_types = FALSE)

if (nrow(df) == 0) {
    message("Table d'enrichissement vide. Figure vide produite.")

    png(outfile, width = 1600, height = 1000, res = 150)
    plot.new()
    title("Analyse d'enrichissement fonctionnel")
    text(0.5, 0.5, "Aucun terme GO enrichi", cex = 1.5)
    dev.off()

    quit(save = "no", status = 0)
}

required_cols <- c("GO_ID", "pvalue", "padj")

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
    stop(
        "Colonnes manquantes : ",
        paste(missing_cols, collapse = ", ")
    )
}

plot_df <- df %>%
    filter(!is.na(pvalue)) %>%
    mutate(
        minus_log10_pvalue = -log10(pvalue),
        label = GO_ID,
        significance = ifelse(padj < 0.05, "padj < 0.05", "padj >= 0.05")
    ) %>%
    arrange(pvalue) %>%
    slice_head(n = top_n)

if (nrow(plot_df) == 0) {
    message("Aucun terme GO avec pvalue valide.")

    png(outfile, width = 1600, height = 1000, res = 150)
    plot.new()
    title("Analyse d'enrichissement fonctionnel")
    text(0.5, 0.5, "Aucun terme GO avec pvalue valide", cex = 1.5)
    dev.off()

    quit(save = "no", status = 0)
}

plot_df <- plot_df %>%
    mutate(label = factor(label, levels = rev(label)))

p <- ggplot(plot_df, aes(x = minus_log10_pvalue, y = label)) +
    geom_col(fill = "darkgreen", color = "black") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    labs(
        title = "Analyse d'enrichissement fonctionnel des gènes associés aux DMR",
        subtitle = "Comparaison : __COMP__ — aphid — __MODTYPE__",
        x = "-log10 P valeur",
        y = NULL
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )

ggsave(
    filename = outfile,
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
)

message("Figure écrite : ", outfile)
'''

        r_code = r_code.replace("__INPUT__", str(input.enrichment))
        r_code = r_code.replace("__OUTPUT__", str(output.plot))
        r_code = r_code.replace("__TOP_N__", str(params.top_n))
        r_code = r_code.replace("__COMP__", str(wildcards.comp))
        r_code = r_code.replace("__MODTYPE__", str(wildcards.modtype))

        with open(r_script, "w") as f:
            f.write(r_code)

        cmd = f"""
        set -euo pipefail
        module load bioconductor 2>/dev/null || true
        Rscript {r_script} > {log[0]} 2>&1
        """

        subprocess.check_call(cmd, shell=True, executable="/bin/bash")