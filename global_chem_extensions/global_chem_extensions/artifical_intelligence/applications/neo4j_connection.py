import pandas as pd
import time
from neo4j import GraphDatabase, basic_auth

class Neo4jConnection:

    def __init__(self, uri, user, pwd):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            print("Failed to create the driver:", e)

    def close(self):
        if self.__driver is not None:
            self.__driver.close()

    def query(self, query, parameters=None, db=None):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        try:
            session = self.__driver.session(database=db) if db is not None else self.__driver.session()
            response = list(session.run(query, parameters))
        except Exception as e:
            print("Query failed:", e)
        finally:
            if session is not None:
                session.close()
        return response

def add_categories(categories):
    # Adds category nodes to the Neo4j graph.
    query = '''
            UNWIND $rows AS row
            MERGE (c:Category {category: row.category})
            RETURN count(*) as total
            '''
    return conn.query(query, parameters = {'rows':categories.to_dict('records')})


def add_names(rows, batch_size=10000):
    # Adds author nodes to the Neo4j graph as a batch job.
    query = '''
            UNWIND $rows AS row
            MERGE (:Name {name: row.name})
            RETURN count(*) as total
            '''
    return insert_data(query, rows, batch_size)


def insert_data(query, rows, batch_size = 10000):
    # Function to handle the updating the Neo4j database in batch mode.

    total = 0
    batch = 0
    start = time.time()
    result = None

    while batch * batch_size < len(rows):

        res = conn.query(query,
                         parameters = {'rows': rows[batch*batch_size:(batch+1)*batch_size].to_dict('records')})
        total += res[0]['total']
        batch += 1
        result = {"total":total,
                  "batches":batch,
                  "time":time.time()-start}
        print(result)

    return result

def add_papers(rows, batch_size=5000):

    query = '''
   UNWIND $rows as row
   MERGE (p:Molecule {id:row.ids}) ON CREATE SET p.smiles = row.smiles
 
   // connect categories
   WITH row, p
   UNWIND row.categories AS category_name
   MATCH (c:Category {category: category_name})
   MERGE (p)-[:IN_CATEGORY]->(c)
 
   // connect authors
   WITH distinct row, p // reduce cardinality
   UNWIND row.names_list AS names
   MATCH (a:Name {name: names})
   MERGE (a)-[:NAMED]->(p)
   RETURN count(distinct p) as total
   '''

    return insert_data(query, rows, batch_size)

if __name__ == '__main__':

    df = pd.read_csv(
        '/Users/sulimansharif/projects/global-chem/global_chem/global_chem_outputs/global_chem.tsv',
        sep='\t',
        names=['names_list', 'smiles', 'node', 'category', 'path']
    )


    paths = df['path'].to_list()

    categories = []
    parent_nodes = []
    ids = []

    global_chem_counters = []
    counter = 1

    for path in paths:

        parent_node  = path.split('.')[-1]
        categories.append(path.split('.')[:-1])
        parent_nodes.append(parent_node)
        global_chem_counters.append(counter)
        counter += 1

    df['nodes'] = parent_nodes
    df['categories'] = categories
    df['ids'] = global_chem_counters

    conn = Neo4jConnection(
        uri="bolt://44.210.87.223:7687",
        user="neo4j",
        pwd="convulsion-endings-bulbs"
    )

    conn.query('CREATE CONSTRAINT molecules IF NOT EXISTS ON (p:Molecule) ASSERT p.id IS UNIQUE')
    conn.query('CREATE CONSTRAINT names IF NOT EXISTS ON (a:Name) ASSERT a.name IS UNIQUE')
    conn.query('CREATE CONSTRAINT categories IF NOT EXISTS ON (c:Category) ASSERT c.category IS UNIQUE')

    categories = pd.DataFrame(df[['categories']])
    categories.rename(columns={'categories':'category'},
                      inplace=True)

    categories = categories.explode('category') \
        .drop_duplicates(subset=['category'])

    names = pd.DataFrame(df[['names_list']])

    names.rename(columns={'names_list':'name'},
                   inplace=True)

    names= names.explode('name').drop_duplicates(subset=['name'])
    add_categories(categories)
    add_names(names)
    add_papers(df)