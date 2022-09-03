/* API to retrieve categories and data from the global-chem repo */

const url = 'https://raw.githubusercontent.com/Sulstice/global-chem/production/global_chem/global_chem_outputs/global_chem.tsv';

const getGlobalChemSmiles = () => {
    return fetch(url, {
        method: 'GET'
    }).then((response) => {
        console.log(response);
        return response.text();
    });
}



export default getGlobalChemSmiles;
