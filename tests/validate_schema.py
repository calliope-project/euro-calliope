import argparse

import jsonschema
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("schema", help="JSON schema in YAML format")
parser.add_argument(
    "--config",
    help=(
        "configuration file to validate. "
        "If not given, schema itself will be validated against JSON schema Draft 2020-12"
    ),
)
parser.add_argument(
    "--level",
    help=(
        "level of the configuration to validate. "
        "If not given, this will validate the entire configuration."
    ),
)
args = parser.parse_args()

with open(args.schema) as f:
    schema = yaml.safe_load(f)

if args.config:
    with open(args.config) as f:
        config = yaml.safe_load(f)
    if args.level:
        config = config[args.level]
    jsonschema.validate(config, schema)
else:
    jsonschema.Draft202012Validator.check_schema(schema)
