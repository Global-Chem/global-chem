import rdflib
import csv
import pandas as pd

if __name__ == '__main__':

    df = pd.read_csv('/Users/sulimansharif/projects/global-chem/global_chem/global_chem.tsv', delimiter='\t', header=None, names=['name', 'smiles', 'node', 'predicate', 'path'])

    # Create the graph object which holds the triples
    graph = rdflib.Graph()

    for i, row in df.iterrows():
        s = rdflib.URIRef(f'#/{row["name"]}')
        p = rdflib.URIRef("#connectsTo")
        o = rdflib.URIRef(f'#/{row["node"]}')
        graph.add((s, p, o))

    for i, row in df.iterrows():

      predicate = row['predicate']

      if str(predicate) == 'nan':
         predicate = 'miscellaenous'

      s = rdflib.URIRef(f'#/{row["node"]}')
      p = rdflib.URIRef("#connectsTo")
      o = rdflib.URIRef(f'#/{predicate}')
      graph.add((s, p, o))

    for i, row in df.iterrows():

      predicate = row['predicate']

      if str(predicate) == 'nan':
        predicate = 'miscellaenous'

      s = rdflib.URIRef(f'#/{predicate}')
      p = rdflib.URIRef("#connectsTo")
      o = rdflib.URIRef(f'#/{"global-chem"}')

      graph.add((s, p, o))

    graph.serialize(destination='graph.ttl', format='application/rdf+xml')
