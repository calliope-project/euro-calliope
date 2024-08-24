"""Creates Calliope top-level model file."""

from pathlib import Path

from eurocalliopelib.template import parametrise_template


def construct_model(
    path_to_template,
    path_to_output,
    modules,
    resolution,
    math_files,
    year,
):
    input_files = sorted([update_path(file, resolution) for file in modules])
    math_files = [update_path(file, resolution) for file in math_files]

    return parametrise_template(
        path_to_template,
        path_to_output,
        input_files=input_files,
        math_files=math_files,
        resolution=resolution,
        year=year,
    )


def update_path(path_string, resolution):
    path = Path(path_string)
    split_point = path.parts.index(resolution)
    return Path(".").joinpath(*path.parts[split_point + 1 :])


if __name__ == "__main__":
    construct_model(
        path_to_template=snakemake.input.template,
        modules=snakemake.input.modules,
        math_files=snakemake.input.math_files,
        resolution=snakemake.wildcards.resolution,
        year=snakemake.params.year,
        path_to_output=snakemake.output[0],
    )
