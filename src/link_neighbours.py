"""Link all those locations that are neighbours."""
import geopandas as gpd
import shapely
import jinja2
import networkx as nx

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
K_EDGE_CONNECTION = 1 # every component of the graph is connect with at least this amount of edges
TEMPLATE = """
links:
    {% for location1, location2 in graph.edges %}
    {{ location1 | replace(".", "-") }},{{ location2 | replace(".", "-") }}:
        techs:
            ac_transmission:
    {% endfor %}
"""


def link_neighbours(path_to_locations, sea_connections, path_to_result):
    """Link all those locations that are neighbours."""
    graph = _create_graph_with_land_connections(path_to_locations)
    if sea_connections:
        graph.add_edges_from([(loc1, loc2) for loc1, loc2 in sea_connections])
    assert nx.is_connected(graph), "There are electrical islands in the network graph."
    links = jinja2.Template(TEMPLATE).render(
        graph=graph
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _create_graph_with_land_connections(path_to_locations):
    regions = gpd.read_file(path_to_locations).to_crs(EPSG_3035_PROJ4).set_index("id")
    graph = nx.Graph()
    graph.add_nodes_from([
        (region_id, {"centroid": region_centroid})
        for region_id, region_centroid in regions.centroid.to_dict().items()
    ])
    neighbours = {
        index: _neighbours(region.geometry, index, regions) for index, region in regions.iterrows()
    }
    graph.add_edges_from([
        (region, other)
        for region, all_others in neighbours.items()
        for other in all_others
    ])
    return graph


def _neighbours(region_geometry, region_index, regions):
    region_geometry = shapely.prepared.prep(region_geometry)
    return [other_index for other_index, other_region in regions.iterrows()
            if (other_index is not region_index) and region_geometry.intersects(other_region.geometry)]


if __name__ == "__main__":
    link_neighbours(
        path_to_locations=snakemake.input.units,
        sea_connections=snakemake.params.sea_connections,
        path_to_result=snakemake.output[0]
    )
