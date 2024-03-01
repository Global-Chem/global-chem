//Render html stuff based on globalchem data

const renderCategories = (data) => {
    let root = document.createElement('ul');
    root.classList.add('list-unstyled', 'ps-0');
    let mainCategory = document.createElement('li');
    mainCategory.classList.add('mb-1');
    let subCategory = document.createElement('div');
    let subCategoryList = document.createElement('ul');
    subCategory.appendChild(subCategoryList);
    const topMenu = document.getElementById('root-ul');
    createMenu(topMenu, data, '', data);

}

const createMenu = (root, menu, category, data) => {
    if (!menu) return;
    let mainCategory = document.createElement('ul');
    mainCategory.classList.add('btn-toggle-nav', 'list-unstyled', 'fw-normal', 'pb-1', 'smaller');
    Object.keys(menu).forEach(submenu => {
        let title = submenu.replaceAll('_', ' ');
        category += `.${submenu}`;
        if (submenu != 'values'){
            let randomValue = Math.floor(Math.random() * 1000);
            let li = document.createElement('li');
            li.classList.add('ms-2');
            let subCategory = document.createElement('div');
            subCategory.classList.add('collapse', 'hide');
            subCategory.id = `${submenu + randomValue}-collapse`;
            let firstLevel = '<button ' +
            'class="btn btn-toggle align-items-center rounded collapsed"' +
            'data-bs-toggle="collapse"' +
            `data-bs-target="#${submenu + randomValue}-collapse"` +
            `aria-expanded="false" >` + title +
            `</button><a href="#${category}" actual="${submenu}" class="badge bg-secondary study-all">Study All</a>`

            if (menu[submenu] && submenu !== 'values'){
                createMenu(subCategory, menu[submenu], category, data);
            }
            if (!menu[submenu].values){
                li.insertAdjacentHTML('afterbegin', firstLevel);
                li.appendChild(subCategory);
            }
            else{
                let a = `<a href="#${category}" class="link-dark rounded">${title}</a>`;
                li.insertAdjacentHTML('afterbegin', a);
            }
            category = category.substr(0, category.lastIndexOf("."));
            mainCategory.appendChild(li);
            if (data.hasOwnProperty(submenu)){
                category = '';
            }
        }
    });
    root.append(mainCategory);
}


export {renderCategories};