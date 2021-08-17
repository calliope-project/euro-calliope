import io
from contextlib import redirect_stdout
from pathlib import Path

import snakemake
import mkdocs
from mkdocs.plugins import BasePlugin
from mkdocs.structure.files import File
import pydot


class DAGPlugin(BasePlugin):
    config_scheme = (
        ('path_to_snakefile', mkdocs.config.config_options.Type(str)),
        ('path_to_src_dir', mkdocs.config.config_options.Type(str)),
        ('path_to_png_relative_to_site', mkdocs.config.config_options.Type(str))
    )

    def on_files(self, files, config, **kwargs):
        """Generate DAG as png and add it to mkdocs files."""
        path_to_snakemake = Path.cwd() / self.config["path_to_snakefile"]
        path_to_png = Path.cwd() / self.config["path_to_src_dir"] / self.config["path_to_png_relative_to_site"]
        path_to_png.parent.mkdir(parents=True, exist_ok=True)
        graph_string = io.StringIO()
        with redirect_stdout(graph_string):
            snakemake.snakemake(path_to_snakemake, printrulegraph=True)
        graph = pydot.graph_from_dot_data(graph_string.getvalue())[0]
        graph.write_png(path_to_png)

        dag_file = File(
            path=self.config["path_to_png_relative_to_site"],
            src_dir=self.config["path_to_src_dir"],
            dest_dir=config["site_dir"],
            use_directory_urls=config["use_directory_urls"]
        )
        files.append(dag_file)
        return files
