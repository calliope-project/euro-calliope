import io
import yaml
from datetime import datetime

import eurocalliopelib


def metadata(config, version, path_to_output):
    metadata = {
        "description":
            "This is the metadata of the build process of "
            "the euro-calliope model in the same directory.",
        "euro-calliope-version": version,
        "euro-calliope-lib-version": eurocalliopelib.__version__,
        "generated-utc": datetime.utcnow(),
        "config": config
    }
    with io.open(path_to_output, 'w', encoding='utf8') as outfile:
        yaml.dump(metadata, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)


if __name__ == "__main__":
    metadata(
        config=snakemake.params.config,
        version=snakemake.params.version,
        path_to_output=snakemake.output[0]
    )
