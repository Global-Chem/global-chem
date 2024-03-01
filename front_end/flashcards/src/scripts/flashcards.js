import { loadGlobalChemValuesFromDisk } from "./disk.js";

let flashCardCounter = {
    index: 0
}

const cardBackTitle = document.getElementById('chem-name');
const cardBackSmiles = document.getElementById('smiles-string');
const categoryTitle = document.getElementById('category-title');


const playFlashcards = (globalChemData) => {
    startGame(globalChemData);
};

const startGame = (globalChemData) => {

    window.addEventListener('hashchange', function () {
        const category = window.location.hash.substring(2);
        let categoryName = window.location.hash;
        categoryName = categoryName.replaceAll('_', ' ').substr(categoryName.lastIndexOf("."), categoryName.length);
        categoryTitle.innerHTML = categoryName.slice(1);

        //Get path
        let smiles = getPropertyPath(globalChemData, category);
        //Randomize flashcards
        smiles = fyShuffle(smiles);
        loadCategory(smiles);
    });
    
}

const loadValuesHelper = (obj, result) => {
    for (const prop in obj) {
        const value = obj[prop];
        if (!value.values) {
            return loadValuesHelper(value, result);
        }
        else {
            result = result.concat(value.values);
        }
    }
    return result;
}

const getPropertyPath = (obj, path) => {
    if (!path){
        if (obj.hasOwnProperty('values')){
            return obj.values;
        }
        //Get all of the nested values arrays
        return loadValuesHelper(obj, []);
    } 
    const properties = path.split('.');
    return getPropertyPath(obj[properties.shift()], properties.join('.'));
}

const fyShuffle = (array) => {
    let currentIndex = array.length,  randomIndex;

    while (currentIndex != 0) {
        randomIndex = Math.floor(Math.random() * currentIndex);
        currentIndex--;

        [array[currentIndex], array[randomIndex]] = [array[randomIndex], array[currentIndex]];
    }

    return array;
}

const loadCategory = (data) => {
    drawSmiles(data[0].smiles);
    changeBack(data[0].name, data[0].smiles);
    const previousButton = document.getElementById('previous');
    const nextButton = document.getElementById('next');
    previousButton.addEventListener('click', () => loadPrevious(data));
    nextButton.addEventListener('click', () => loadNext(data));
}

const loadNext = (data) => {
    if (flashCardCounter.index === data.length - 1){
        flashCardCounter.index = 0;
    }
    else{
        flashCardCounter.index += 1;
    }
    let index = flashCardCounter.index;
    changeBack(data[index].name, data[index].smiles);
    drawSmiles(data[index].smiles);
}

const changeBack = (name, smiles) => {
    cardBackTitle.innerText = name;
    cardBackSmiles.innerText = smiles;
}

const loadPrevious = (data) => {
    if (flashCardCounter.index === 0){
        flashCardCounter.index = data.length - 1;
    }
    else{
        flashCardCounter.index -= 1;
    }
    let index = flashCardCounter.index;
    changeBack(data[index].name, data[index].smiles);
    drawSmiles(data[index].smiles);
}

const drawSmiles = (smile) => {
    let options = {width: 200, height: 200};

    // Initialize the drawer to draw to canvas
    let smilesDrawer = new SmilesDrawer.Drawer(options);
    SmilesDrawer.parse(smile, function(tree) {
        // Draw to the canvas
        smilesDrawer.draw(tree, "example-canvas", "light", false);
        // Alternatively, draw to SVG:
        // svgDrawer.draw(tree, 'output-svg', 'dark', false);
    });
    // Alternatively, initialize the SVG drawer:
    // let svgDrawer = new SmilesDrawer.SvgDrawer(options);

}

export {playFlashcards};