import yaml
from pathlib import Path

import jsonschema2md
import mkdocs
from mkdocs.plugins import BasePlugin
from mkdocs.structure.files import File

DESCRIPTION = """
Below you can find a complete enumeration of all configuration parameters of Euro-Calliope's
workflow including short descriptions and datatypes. To learn how to change the parameter
values for your model builds, head over to the Customisation section of the workflow documentation.
"""


class SchemaPlugin(BasePlugin):
    config_scheme = (
        ('path_to_schema', mkdocs.config.config_options.Type(str)),
        ('path_to_src_dir', mkdocs.config.config_options.Type(str)),
        ('path_to_md_relative_to_site', mkdocs.config.config_options.Type(str))
    )

    def on_serve(self, server, config, builder):
        path_to_schema = Path.cwd() / self.config["path_to_schema"]
        server.watch(path_to_schema)
        return server

    def on_files(self, files, config, **kwargs):
        """Generate overview over configuration schema and add to mkdoc's files."""
        path_to_schema = Path.cwd() / self.config["path_to_schema"]
        path_to_md = Path.cwd() / self.config["path_to_src_dir"] / self.config["path_to_md_relative_to_site"]
        path_to_md.parent.mkdir(parents=True, exist_ok=True)

        parser = jsonschema2md.Parser()
        with path_to_schema.open("r") as f_schema:
            schema = yaml.safe_load(f_schema)
        schema = SchemaPlugin.move_quality_control(schema)
        with path_to_md.open("w") as f_md:
            lines = parser.parse_schema(schema)
            lines = SchemaPlugin.customise_markdown(lines)
            f_md.writelines(lines)

        schema_md_file = File(
            path=self.config["path_to_md_relative_to_site"],
            src_dir=self.config["path_to_src_dir"],
            dest_dir=config["site_dir"],
            use_directory_urls=config["use_directory_urls"]
        )
        files.append(schema_md_file)
        return files

    @staticmethod
    def move_quality_control(schema):
        # Move quality-control to the end as it causes issues
        qc = schema['properties']['quality-control']
        del schema['properties']['quality-control']
        schema['properties']['quality-control'] = qc
        return schema

    @staticmethod
    def customise_markdown(lines):
        # 1. Move from 2 spaces to 4 spaces intend as the 2 spaces don't work properly with mkdocs
        lines = ['    '.join(line.split('  ')) for line in lines]
        # 2. Change headline
        assert lines[0] == '# JSON Schema\n\n'
        lines[0] = "# Configuration parameters of Euro-Calliope's workflow\n\n"
        # 3. Remove main description and subheadline
        assert lines[2] == '## Properties\n\n'
        del lines[2]
        #del lines[1]
        return lines
