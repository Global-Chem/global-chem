//Saves global chem data to local storage to limit fetch (api) requests


const saveGlobalChemToDisk = (parsedData) => {
    const jsonedData = JSON.stringify(parsedData);
    try {
        window.localStorage.setItem('globalChemData', jsonedData);
        window.localStorage.setItem('lastUpdated', new Date().getTime());
        return true;
    } catch (e){
        console.log(e);
        return e;
    }
}

const loadGlobalChemFromDisk = () => {
    try {
        let globalChemData = window.localStorage.getItem('globalChemData');
        //Turn it back into JSON object
        let parsedData = JSON.parse(globalChemData);
        return parsedData;
    } catch (e){
        console.log(e);
        return e;
    }
}

const loadGlobalChemValuesFromDisk = () => {
    try {
        let globalChemValues = window.localStorage.getItem('globalChemValues');
        //Turn it back into JSON object
        let parsedData = JSON.parse(globalChemValues);
        return parsedData;
    } catch (e){
        console.log(e);
        return e;
    }
}

const loadLastUpdatedFromDisk = () => {
    try {
        let lastUpdated = window.localStorage.getItem('lastUpdated');
        return lastUpdated;
    } catch (e){
        console.log(e);
        return e;
    }
}

export {saveGlobalChemToDisk, loadGlobalChemFromDisk, loadLastUpdatedFromDisk, loadGlobalChemValuesFromDisk}