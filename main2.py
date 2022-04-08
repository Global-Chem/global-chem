
from global_chem.global_chem import GlobalChem

from urllib.parse import quote

import requests

if __name__ == '__main__':

    url = 'https://chemwriter.com/smiles/?smiles=%s' % quote('C1CCC1')
    response = requests.get(url)
    print (response.content.decode('utf-8'))