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
        ('path_to_png', mkdocs.config.config_options.Type(str)),
        ('path_to_cache', mkdocs.config.config_options.Type(str))
    )

    def on_pre_build(self, config, **kwargs):
        """Generate DAG as png."""
        path_to_snakemake = Path.cwd() / self.config["path_to_snakefile"]
        path_to_png = Path.cwd() / self.config["path_to_png"]
        path_to_cache = Path.cwd() / self.config["path_to_cache"]
        path_to_png.parent.mkdir(parents=True, exist_ok=True)
        path_to_cache.parent.mkdir(parents=True, exist_ok=True)
        graph_string = io.StringIO()
        with redirect_stdout(graph_string):
            snakemake.snakemake(path_to_snakemake, printrulegraph=True)
        graph = pydot.graph_from_dot_data(graph_string.getvalue())[0]
        try:
            previous_graph = pydot.graph_from_dot_file(path_to_cache)[0]
        except (FileNotFoundError, TypeError):
            previous_graph = pydot.Dot()
        # Only generate PNG if it does not exist or is outdated.
        # This is important for mkdoc's livereload feature, which will constantly
        # reload otherwise.
        if path_to_png.exists() and (graph.to_string() != previous_graph.to_string()):
            graph.write_png(path_to_png)
            graph.write_raw(path_to_cache)
