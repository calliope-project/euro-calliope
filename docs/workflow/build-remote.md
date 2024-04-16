# Build on a cluster

You may want to build the model on a cluster.
While you can build Euro-Calliope on [any cluster that is supported by Snakemake](https://snakemake.readthedocs.io/en/v{{ snakemake_version }}/tutorial/additional_features.html#cluster-or-cloud-execution), our default configuration is targeted at, and tested on, ETH's Euler cluster.

## Build

To build the model on Euler, use the following command:

```bash
snakemake --profile profiles/euler
```

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/v{{ snakemake_version }}/tutorial/additional_features.html#cluster-or-cloud-execution) and take `profiles/euler` as a starting point.

## Work local, build on remote

If you are like us, you may want to work locally (to change configuration parameters, add modules etc), but execute remotely on the cluster.
We support this workflow through three Snakemake rules: `send`, `receive`, and `clean_cluster_results`.
It works like the following.

First, start local and make sure the `cluster-sync` configuration parameters fit your environment.
Next, run `snakemake send` to send the entire repository to your cluster.
On the cluster, execute the workflow with Snakemake ([see above](./build-remote.md#build-on-a-cluster)).
After the workflow has finished, download results by locally running `snakemake receive`.
By default, this will download results into `build/cluster`.

This workflow works iteratively too.
After analysing your cluster results locally, you may want to make changes locally, send these changes to the cluster (`snakemake send`), rerun on the cluster, and download updated results (`snakemake receive`).

To remove cluster results on your local machine, run `snakemake clean_cluster_results`.

## Be notified of build successes or fails

 As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully.
 Notifications are sent by email.
 To activate notifications, add the email address of the recipient to the configuration key `email`.
 You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

```bash
snakemake --profile profiles/euler --config email=<your-email>
```
