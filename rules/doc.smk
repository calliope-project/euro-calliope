"""Rules for generating documentation artefacts."""

include: "./shapes.smk"
configfile: "./config/default.yaml"


rule spatial_scope_and_resolutions:
    message: "Plot spatial scope and resolution for documentation."
    input:
        script = "scripts/vis/spatial_scope_and_resolutions.py",
        regional_units = "build/data/regional/units.geojson",
        national_units = "build/data/national/units.geojson"
    params:
        dpi = 300
    conda: "../envs/vis.yaml"
    output: "docs/img/spatial-scope-and-resolutions.png"
    script: "../scripts/vis/spatial_scope_and_resolutions.py"
