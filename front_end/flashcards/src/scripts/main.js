import getGlobalChemSmiles from './api.js';
import parseTSV from './parser.js';
import {
    saveGlobalChemToDisk, 
    loadGlobalChemFromDisk, 
    loadLastUpdatedFromDisk
} from './disk.js'
import {renderCategories} from './render.js';
import { playFlashcards } from './flashcards.js';

//Check if localstorage contains globalchem data already

const main = async () => {
    const globalChemData = await loadGlobalChemData();
    renderCategories(globalChemData);
    playFlashcards(globalChemData);
}

const loadGlobalChemData = async () => {
    const currentDate = new Date().getTime();
    const lastUpdated = loadLastUpdatedFromDisk();
    const millisecondsCheck = 3600000; //An hour in milliseconds
    
    //New user, nothing stored in local storage or it has been an hour or more since last update
    if (!lastUpdated || currentDate - lastUpdated > millisecondsCheck){
        const data = await getGlobalChemSmiles();
        const parsedData = parseTSV(data);
        saveGlobalChemToDisk(parsedData);
        return parsedData;
    }
    
    //Retrieve from local storage
    else{
        const globalChemData = loadGlobalChemFromDisk();
        return globalChemData;
    }
    
}

main();

