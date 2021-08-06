import io
from contextlib import redirect_stdout
from pathlib import Path

import snakemake
import mkdocs
from mkdocs.plugins import BasePlugin
import pydot


class DAGPlugin(BasePlugin):
    config_scheme = (
        ('path_to_snakefile', mkdocs.config.config_options.Type(str)),
        ('path_to_png', mkdocs.config.config_options.Type(str))
    )

    def on_pre_build(self, config):
        path_to_snakemake = Path.cwd() / self.config["path_to_snakefile"]
        path_to_png = Path.cwd() / self.config["path_to_png"]
        path_to_png.parent.mkdir(parents=True, exist_ok=True)
        graph_string = io.StringIO()
        with redirect_stdout(graph_string):
            snakemake.snakemake(path_to_snakemake, printrulegraph=True)
        graph = pydot.graph_from_dot_data(graph_string.getvalue())
        graph[0].write_png(path_to_png)
