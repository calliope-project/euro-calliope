localrules: eurostat_data_tsv

rule eurostat_data_tsv:
    message: "Get {wildcards.dataset} from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"]["eurostat-base-url"].format(dataset=wildcards.dataset)
    output: protected("data/automatic/eurostat-{dataset}.tsv.gz")
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat data"
    input:
        src = script_dir + "eurostat/annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat-nrg_bal_c.tsv.gz",
    output: "build/data/eurostat/annual-energy-balances.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"
