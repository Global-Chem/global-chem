/* TSV parser for category and smiles data */

//Input: TSV as a string 
//Output: Object with correct parent and nested categories with smiles values
const parseTSV = (tsv) => {
    let globalChemOrganizedData = {};
    //Split into lines
    const splitLines = (tsv.split(/\n/)).slice(0, -1); //Remove empty character
    for (let i = 0; i < splitLines.length; i++){
        //Split each line by tabs
        const splitTabs = splitLines[i].split(/\t/);
        //Third index is the directory (category) path for each line
        const category = splitTabs[4].slice(12); //Minus global_chem start path
        const name = splitTabs[0];
        const smiles = splitTabs[1];

        addData(globalChemOrganizedData, category, {name, smiles});

    }
    console.log(globalChemOrganizedData)
    return globalChemOrganizedData;
}

const addData = (obj, path, val) => { 
    //Splits path into individual props
    const keys = path.split('.');
    const lastKey = 'values';
    const lastObj = keys.reduce((obj, key) => {
        return obj[key] = obj[key] || {}
    }, obj);
    lastObj[lastKey] === undefined ? lastObj[lastKey] = [val] : lastObj[lastKey].push(val);
    return lastObj;
};



export default parseTSV;