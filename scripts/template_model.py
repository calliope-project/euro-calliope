"""Creates Calliope top-level model file."""
from eurocalliopelib.parametrise_template import parametrise_template


def construct_model(path_to_template, path_to_output, config_files, resolution, model_year):
    config_files = sorted([
        file.replace(f"build/models/{resolution}", ".") for file in config_files
    ])

    return parametrise_template(
        path_to_template, path_to_output,
        config_files=config_files,
        resolution=resolution,
        model_year=model_year
    )


if __name__ == "__main__":
    construct_model(
        path_to_template=snakemake.input.template,
        config_files=snakemake.input.config_files,
        resolution=snakemake.wildcards.resolution,
        model_year=snakemake.params.model_year,
        path_to_output=snakemake.output[0],
        )
