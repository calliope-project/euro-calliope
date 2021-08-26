import mkdocs
from mkdocs.plugins import BasePlugin
from mkdocs.structure.files import File


class ChangelogPlugin(BasePlugin):
    config_scheme = (
        ('path_to_changelog', mkdocs.config.config_options.Type(str)),
        ('path_to_src_dir', mkdocs.config.config_options.Type(str))
    )

    def on_files(self, files, config, **kwargs):
        """Link top-level changelog to mkdoc's files."""
        changelog_md_file = File(
            path=self.config["path_to_changelog"],
            src_dir=self.config["path_to_src_dir"],
            dest_dir=config["site_dir"],
            use_directory_urls=config["use_directory_urls"]
        )
        files.append(changelog_md_file)
        return files
