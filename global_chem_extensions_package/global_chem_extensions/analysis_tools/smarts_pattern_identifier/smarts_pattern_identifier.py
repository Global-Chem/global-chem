#!/usr/bin/env python3
#
# GlobalChemExtensions - SmartsPatternIdentifier
#
# ----------------------------------------------

# Imports
# -------

import sys
import json

# Flask Imports
# -------------

from flask import Flask
from flask import request
from flask import render_template_string

# RDkit Imports
# -------------

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions

# Drawing Options
# ---------------

DrawingOptions.atomLabelFontSize = 55
DrawingOptions.dotsPerAngstrom = 100
DrawingOptions.bondLineWidth = 2.0

# MiniFrag Database
# -----------------

# TODO: Move this into GlobalChem

minifrag_database = [
    'Brc1cn[nH]c1',
    'C1COCCN1',
    'c1cn[nH]c1',
    'OCC1CC1',
    'C1CCOC1',
    'Nc1ccon1',
    'O=C(O)[C@@H]1CCN1',
    'N#Cc1ccon1',
    'OC1CNC1',
    'O=P(O)(O)O',
    'O=C1CCON1',
    'N#Cc1c[nH]cn1',
    'CCC(=O)NO',
    'CN1CC[C@H](O)C1',
    'Nc1ccccn1',
    'N#CC1CC1',
    'NC(=O)C1(N)CC1',
    'O=S1(=O)CC(O)C1',
    'Nc1cc[nH]n1',
    'NC1CCCC1',
    'CC(C)/C(N)=N/O',
    'O=S1(=O)CCC1',
    'Oc1ncccn1',
    'NC1=NCCCN1',
    'O=C1CCC(=O)N1',
    'OCc1cc[nH]n1',
    'Nc1ncc[nH]1',
    'c1c[nH]cn1',
    'Nc1csc(N)n1',
    'N[C@@H]1CCNC1=O',
    'O=C1CNCCN1',
    'CC(C)NC(N)=O',
    'Nc1cccnn1',
    'OC1CCCCC1',
    'Nc1nc[nH]n1',
    'NCc1cocn1',
    'C[C@H](O)C(=O)O',
    'c1ccncc1',
    'O=c1cccn[nH]1',
    'Oc1ccnnc1',
    'c1c[nH]nn1',
    'CNC(=O)[C@H](C)N',
    'OB(O)C1CCC1',
    'N#CCC(N)=O',
    'Oc1ccccc1',
    'C1CCNC1',
    'NCC(=O)O',
    'CC[C@H](N)C(N)=O',
    'OC[C@@H]1CCCN1',
    'C#C[C@H](C)N',
    'O=C1CCNN1',
    'Cc1nnn[nH]1',
    'O=C1CCCCN1',
    'CCC(=O)O',
    'c1nnco1',
    'CC(=O)NC(N)=O',
    'O=C(O)C1CCC1',
    'O=C1NCCN1',
    'CC(=O)C(=O)O',
    'Nc1ccncc1',
    'C1CNCCN1',
    'O=C(O)c1ccsc1',
    'O=C1CNC(=O)N1',
    'NCc1ccsc1',
    'Nc1ncccn1',
    'Brc1cccnc1',
    'O[C@H]1COC[C@H]1O',
    'C[C@H](N)CC(N)=O',
    'Cn1cnc(N)c1',
    'Nc1cc(O)n[nH]1',
    'C1CNCCOC1',
    'NCC1CC1',
    'O=C1CCCN1',
    'Oc1ccccn1',
    'O=S1(=O)CCN1',
    'CC(C)CCCN',
    'Nc1cncnc1',
    'CCC(=N)N',
    'CCOC(N)=O',
    'NC(=O)C1CCC1',
    'NC[C@@H]1CCCO1'
]

class SmartsPatternIdentifier(object):

    '''

    Flask Application to run SMARTS Patterns.

    '''

    def __init__(self):

        self.html = self.create_html_template()
        self.app = self.create_app()


    def create_html_template(self):

        '''

        Create the HTML Template

        '''

        import textwrap

        html = '''\
        
        <html>
        <head>
          <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta2/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BmbxuPwQa2lc/FVzBcNJ7UAyJxM6wuqIj61tLrc4wSX0szH/Ev+nYRRuWlolflfl" crossorigin="anonymous">
          <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.21/css/jquery.dataTables.css">
          <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js"   integrity="sha256-4+XzXVhsDmqanXGHaHvgh1gMQKX40OUvDEBTu8JcmNs="   crossorigin="anonymous"></script>
          <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.21/js/jquery.dataTables.js"></script>
        </head>
        <body>
          <div class="container">
        
            <h1>SMARTS Validator</h1>
            <p class="lead">Test SMARTS patterns using set of molecules.</p>
        
            <form method="POST" enctype="multipart/form-data">
              <div class="mb-3">
                <label for="smarts" class="form-label">SMARTS</label>
                <textarea class="form-control" id="smarts" name="smarts" rows="10">
                    {
                        "halogen": "[F,Cl,Br,I]",
                        "acetic anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
                        "hydroxyl": "[OX2H]",
                        "phenol": "[OX2H][cX3]:[c]",
                        "aldehyde": "[CX3H1](=O)[#6]",
                        "amide": "[NX3][CX3](=[OX1])[#6]",
                        "ketone": "[#6][CX3](=O)[#6]",
                        "ester": "[#6]-C(=O)O-[#6]",
                        "ether": "[OD2]([#6])[#6]",
                        "amine": "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",
                        "enamine": "[NX3][CX3]=[CX3]",
                        "vinylic alkene": "[$([CX3]=[CX3])]",
                        "alkyne": "[$([CX2]#C)]",
                        "allenic alkene": "[$([CX2](=C)=C)]",
                        "sulfide": "[#16X2H0]",
                        "sulfonamide": "[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]",
                        "sulfoxide": "[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]"
                    }
                </textarea>
              </div>
              <div class="mb-3">
                <label class="form-label">Test Molecules</label>
                <div class="form-check">
                  <input class="form-check-input" type="radio" name="testfile" id="testfile" value="astex_mini_frag.sdf" checked>
                  <label class="form-check-label" for="testfile">
                    Astex Mini Frag
                  </label>
                </div>
              </div>
              <button type="submit" class="btn btn-primary">Submit</button>
            </form>
        
        {% if data %}
            <div>
              <table class="table" id="table">
                <thead>
                  <tr>
                    <th scope="col">Molecule</th>
                    {% for p in patts %}
                    <th scope="col">{{ p }}</th>
                    {% endfor %}
                  </tr>
                </thead>
                <tbody>
                {% for entry in data %}
                  <tr>
                    <td>{{ entry.mol|safe }}</td>
                    {% for p in patts %}
                    <td>{{ entry[p]|safe }}</td>
                    {% endfor %}
                  </tr>
                {% endfor %}
                </tbody>
              </table>
            </div>
        {% endif %}
          </div>
        }
          <script>
            $(document).ready( function () {
              $('#table').DataTable({
                paging: false
              });
            } );
          </script>
        </body>
        </html>
        '''

        return textwrap.dedent(html)

    def create_app(self):

        '''

        Create the application

        '''

        app = Flask(__name__)

        @app.route('/', methods=['POST', 'GET'])
        def create_smarts_data():

            if request.method == 'POST':

                smarts = json.loads(request.form['smarts'])
                patterns = {key: Chem.MolFromSmarts(value) for key, value in smarts.items()}

                mols = [ Chem.MolFromSmiles(m) for m in minifrag_database ]

                data = []

                for m in mols:

                    tm = Chem.RemoveHs(m)
                    AllChem.Compute2DCoords(tm)
                    drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)
                    drawer.DrawMolecule(tm)
                    drawer.FinishDrawing()
                    svg = drawer.GetDrawingText()
                    entry = {'mol': svg}

                    for name, p in patterns.items():

                        if tm.HasSubstructMatch(p):

                            tm = Chem.RemoveHs(m)
                            AllChem.Compute2DCoords(tm)
                            drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)
                            hs = tm.GetSubstructMatches(p)
                            drawer.DrawMolecule(tm, highlightAtoms=sum(hs, ()))
                            drawer.FinishDrawing()
                            svg = drawer.GetDrawingText()
                            entry[name] = svg

                        else:
                            entry[name] = ''
                    data.append(entry)

                return render_template_string(self.html, data=data, patts=patterns)

            return render_template_string(self.html)

        return app

    def launch_app(self, host='0.0.0.0', port=5000, debug=True):

        '''

        Launch the application

        Arguments:
            host (String): host of where the user would like to host the server
            port (Int): port number the user would like to put the application
            debug (Bool): debugger flag for the flask app.
        '''


        self.app.run(host=host, port=port, debug=debug)
