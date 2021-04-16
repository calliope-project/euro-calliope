import sys
import jsonschema
import yaml

if __name__ == "__main__":
    configfile = sys.argv[1]
    with open(configfile, "r") as f:
        config = yaml.safe_load(f)

    if len(sys.argv) > 2:
        schemafile = sys.argv[2]
        with open(schemafile, "r") as f:
            schema = yaml.safe_load(f)
        jsonschema.validate(config, schema)
    else:
        jsonschema.Draft7Validator.check_schema(config)
