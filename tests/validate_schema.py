import argparse
import jsonschema
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("schema", help="JSON schema in YAML format")
parser.add_argument(
    "--config",
    help=(
        "configuration file to validate. "
        "If not given, schema itself will be validated against JSON schema Draft 7"
    )
)
args = parser.parse_args()

with open(args.schema, "r") as f:
    schema = yaml.safe_load(f)

if args.config:
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    jsonschema.validate(config, schema)
else:
    jsonschema.Draft7Validator.check_schema(schema)
