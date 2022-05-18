import os
import glob
import yaml
from pathlib import Path

import mkdocs
from mkdocs.plugins import BasePlugin
from mkdocs.structure.files import File

from eurocalliopelib.template import parametrise_template


class TechModulePlugin(BasePlugin):
    config_scheme = (
        ('path_to_md_template_to_fill', mkdocs.config.config_options.Type(str)),
        ('path_to_md_relative_to_site', mkdocs.config.config_options.Type(str)),
        ('path_to_tech_module_contents', mkdocs.config.config_options.Type(str)),
        ('path_to_src_dir', mkdocs.config.config_options.Type(str))

    )

    def on_serve(self, server, config, builder):
        path_to_tech_module_contents = Path.cwd() / self.config["path_to_tech_module_contents"]
        server.watch(path_to_tech_module_contents)
        return server

    def on_files(self, files, config, **kwargs):
        """Generate overview over configuration schema and add to mkdoc's files."""
        path_to_template_md = Path.cwd() / self.config["path_to_md_template_to_fill"]
        path_to_tech_module_contents = Path.cwd() / self.config["path_to_tech_module_contents"]
        path_to_md = Path.cwd() / self.config["path_to_src_dir"] / self.config["path_to_md_relative_to_site"]
        path_to_md.parent.mkdir(parents=True, exist_ok=True)

        with open(path_to_tech_module_contents) as f:
            tech_modules = yaml.safe_load(f)

        parametrise_template(path_to_template_md, path_to_md, modules=tech_modules)

        final_md_file = File(
            path=self.config["path_to_md_relative_to_site"],
            src_dir=self.config["path_to_src_dir"],
            dest_dir=config["site_dir"],
            use_directory_urls=config["use_directory_urls"]
        )
        files.append(final_md_file)
        return files
