import sys
from global_chem import g

gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)
network = gc.network
networkx_graph = cheminformatics.convert_to_networkx(network)
print (networkx_graph.nodes.data())


import networkx as nx

from vispy import app, scene
from vispy.visuals.graphs import layouts


canvas = scene.SceneCanvas(title='Simple NetworkX Graph', size=(600, 600),
                           bgcolor='white', show=True)
view = canvas.central_widget.add_view('panzoom')

graph = nx.adjacency_matrix(
  nx.fast_gnp_random_graph(500, 0.005, directed=True))
layout = layouts.get_layout('force_directed', iterations=100)

visual = scene.visuals.Graph(
  graph, layout=layout, line_color='black', arrow_type="stealth",
  arrow_size=30, node_symbol="disc", node_size=20,
  face_color=(1, 0, 0, 0.2), border_width=0.0, animate=True, directed=False,
  parent=view.scene)


@canvas.events.draw.connect
def on_draw(event):
  if not visual.animate_layout():
    canvas.update()

if __name__ == '__main__':
  if sys.flags.interactive != 1:
    app.run()
