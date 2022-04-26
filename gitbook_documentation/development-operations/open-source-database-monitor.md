# Open Source Database Monitor

With the rise of cheminformatics, so comes the rise of open source databases. Unfortunately, open source databases are spread amongst a variety of URLS and can be hard to keep track in an efficient fashion. We create the `Uptime` bot

{% embed url="https://chemistrydb.com" %}

```
 Zinc 15                                                  OpenFDA                                                      
 Zinc 20                                                  Metabolites Biological Role                                  
 PubChem                                                  MetaboAnalyst                                                
 NIST Chemistry Webhook                                   Adverse Drug Reaction Classification System                  
 Chem Exper                                               Metabolism and Transport Database                            
 NMR Shift Database                                       Ecology Toxicity                                             
 Drug Bank                                                Human and Environment Risk Assessment                        
 Binding Database                                         International Toxicity Information for Risk Assesments       
 Spectral Database for Organic Compounds                  Japan Exisiting Database                                     
 Sider                                                    National Pesticide Center                                    
 ChemSpider                                               Pesticide Info                                               
 Stitch                                                   Kyoto Encyclopedia of Genes and Genomes                      
 CardPred                                                 Hetereocycles                                                
 Comparative Toxicogenomics Database                      Chemical Resolver                                            
 AMED Cardiotoxicity Database                             LookChem                                                     
 Tox21                                                    Lipid Maps                                                   
 Drug Safety Analysis System                             
```

**Imports**

```
from global_chem_extensions import GlobalChemExtensions

do = GlobalChemExtensions().development_operations()
```

{% tabs %}
{% tab title="Code" %}
```
successes, failures = do.check_status_on_open_source_databases()
print (successes)
print (failures)
```
{% endtab %}

{% tab title="Output" %}
```
{' Zinc 15': 'Up', ' Zinc 20': 'Up', ' PubChem': 'Up', ' NIST Chemistry Webhook': 'Up', ' Chem Exper': 'Up', ' NMR Shift Database': 'Up', ' Binding Database': 'Up', ' Spectral Database for Organic Compounds': 'Up', ' Sider': 'Up', ' ChemSpider': 'Up', ' Stitch': 'Up', ' CardPred': 'Up', ' Comparative Toxicogenomics Database': 'Up', ' AMED Cardiotoxicity Database': 'Up', ' Tox21': 'Up', ' Drug Safety Analysis System': 'Up', ' OpenFDA': 'Up', ' Metabolites Biological Role': 'Up', ' MetaboAnalyst': 'Up', ' Adverse Drug Reaction Classification System': 'Up', ' Metabolism and Transport Database': 'Up', ' Ecology Toxicity ': 'Up', ' Japan Exisiting Database': 'Up', ' National Pesticide Center': 'Up', ' Pesticide Info': 'Up', ' Kyoto Encyclopedia of Genes and Genomes': 'Up', ' Hetereocycles': 'Up', ' Chemical Resolver': 'Up', ' LookChem': 'Up', ' Lipid Maps': 'Up'}
{' Drug Bank': 'Down', ' Human and Environment Risk Assessment ': 'Down', ' International Toxicity Information for Risk Assesments': 'Down'}
```
{% endtab %}
{% endtabs %}

**Github Repository**

![](<../.gitbook/assets/Screen Shot 2022-02-24 at 10.26.09 AM.png>)

The bot serves as a heart beat checker pinging URLS for databases that are most often used by cheminformaticians.&#x20;

![](<../.gitbook/assets/Screen Shot 2022-02-24 at 10.28.37 AM.png>)

Within the GlobalChemExtensions component we would like a persistent monitoring system tied to the same URLS the bot is also pinging.&#x20;

```
for i in contents:
    if 'name' in i and \
        '#' not in i and 'chemistrydb.com' not in i \
        and 'Cheminformatic Database Statuses' not in i:
        i = i.strip('\n').split('name:')[1]
        names.append(i)

    if 'url' in i:
        i = i.strip('\n').split('url:')[1]
        urls.append(i)

return urls, names
```

The code then checks the status code to see if it still alive. If not it gets reported as a failure. This is useful for folk relying on these databases for whatever processing they are doing.&#x20;

```
response = urllib.request.urlopen(url)
status_code = response.getcode()

if '200' == status_code:
    successes[name] = 'Up'
 except:
    failures[name] = 'Down'
```

