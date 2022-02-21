function getGraphDataSets() {

    const loadMiserables = function(Graph) {
        Graph
            .cooldownTicks(400)
            .nodeLabel('id')
            .nodeAutoColorBy('group')
            .forceEngine('ngraph')
            .jsonUrl('front_end/datasets/chemical_graph_network.json');
    };
    loadMiserables.description = "<em>GlobalChem Graph Network</em>";

    //

    const loadBlocks = function(Graph) {
        fetch('front_end/datasets/blocks.json').then(r => r.json()).then(data => {
            data.nodes.forEach(node => { node.name = `${node.user?node.user+': ':''}${node.description || node.id}` });

        Graph
            .cooldownTicks(300)
            .cooldownTime(10000000)
            .nodeAutoColorBy('user')
            .forceEngine('ngraph')
            .graphData(data);
    });
    };
    loadBlocks.description = "<em>Blocks</em> data (<a href='https://bl.ocks.org/mbostock/afecf1ce04644ad9036ca146d2084895'>afecf1ce04644ad9036ca146d2084895</a>)";

    //

    const loadD3Dependencies = function(Graph) {
        fetch('.d3.csv').then(r => r.text()).then(d3.csvParse).then(data => {
            const nodes = [], links = [];
        data.forEach(({ size, path }) => {
            const levels = path.split('/'),
            module = levels.length > 1 ? levels[1] : null,
            leaf = levels.pop(),
            parent = levels.join('/');

        nodes.push({
            path,
            leaf,
            module,
            size: +size || 1
        });

        if (parent) {
            links.push({ source: parent, target: path});
        }
    });

        Graph
            .cooldownTicks(300)
            .nodeRelSize(0.5)
            .nodeId('path')
            .nodeVal('size')
            .nodeLabel('path')
            .nodeAutoColorBy('module')
            .forceEngine('ngraph')
            .graphData({ nodes: nodes, links: links });
    });
    };
    loadD3Dependencies.description = "<em>D3 dependencies</em> data (<a href='https://bl.ocks.org/mbostock/9a8124ccde3a4e9625bc413b48f14b30'>9a8124ccde3a4e9625bc413b48f14b30</a>)";

    const tunnel = function(Graph) {

        const perimeter = 12, length = 30;

        const getId = (col, row) => `${col},${row}`;

        let nodes = [], links = [];
        for (let colIdx=0; colIdx<perimeter; colIdx++) {
            for (let rowIdx=0; rowIdx<length; rowIdx++) {
                const id = getId(colIdx, rowIdx);
                nodes.push({id});

                // Link vertically
                if (rowIdx>0) {
                    links.push({ source: getId(colIdx, rowIdx-1), target: id });
                }

                // Link horizontally
                links.push({ source: getId((colIdx || perimeter) - 1, rowIdx), target: id });
            }
        }

        Graph
            .cooldownTicks(300)
            .forceEngine('ngraph')
            .graphData({ nodes: nodes, links: links });
    };
    tunnel.description = "fabric data for a cylindrical tunnel shape";

    //

    return [loadMiserables, loadBlocks, loadD3Dependencies, tunnel];
}
