# SMARTS Pattern Identifier

We want to validate our SMARTS strings to verify that it identifies exactly the match of functional group we would like when presented with different scenarios. The MiniFragDatabase is a good example of testing your SMARTS string against that database.

**Imports**

```
from global_chem_extensions import GlobalChemExtensions
cheminformatics = GlobalhemExtensions().cheminformatics()
```

**Load the Application**

```

spi = cheminformatics.smarts_pattern_identifier(
    host='0.0.0.0',
    port=5000,
    debugger=True
)

```

This will launch the app with a fixed amount of SMARTS strings for the user to get familiar with the application.&#x20;

If we look at one example that's built like the hydroxyl group:

```
"hydroxyl": "[OX2H]",
```

We can test to see what environments this SMARTs matches in a MiniFrag Database set.&#x20;

![](<../.gitbook/assets/Screen Shot 2022-03-09 at 10.39.55 AM.png>)

From here we can see that the `OH` hydroxyl group is matching a carboxylic acid as well as a regular hydroxyl group. Which means that you would probably need to adjust the smarts string for a match like that. It would need to be more acute.&#x20;

