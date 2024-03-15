import mkdocs
from mkdocs.plugins import BasePlugin
from mkdocs.structure.files import File


class AddFilePlugin(BasePlugin):
    config_scheme = (
        ("path_to_file", mkdocs.config.config_options.Type(str)),
        ("path_to_src_dir", mkdocs.config.config_options.Type(str)),
    )

    def on_files(self, files, config, **kwargs):
        """Link top-level files to mkdocs files."""
        linked_md_file = File(
            path=self.config["path_to_file"],
            src_dir=self.config["path_to_src_dir"],
            dest_dir=config["site_dir"],
            use_directory_urls=config["use_directory_urls"],
        )
        files.append(linked_md_file)
        return files
