function tissueFormSearch() {
    // Declare variables 
    var input, filter, table, tr, td, i;
    input = document.getElementById("tissueFormNameInput");
    filter = input.value.toUpperCase();
    table = document.getElementById("tissueFormTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[0];
    if (td) {
        if (td.innerHTML.toUpperCase().indexOf(filter) > -1) {
        tr[i].style.display = "";
        } else {
        tr[i].style.display = "none";
        }
    } 
    }
}

function tfFormSearch() {
    // Declare variables 
    var input, filter, table, tr, td, i;
    input = document.getElementById("tfFormNameInput");
    filter = input.value.toUpperCase();
    table = document.getElementById("tfFormTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[0];
    if (td) {
        if (td.innerHTML.toUpperCase().indexOf(filter) > -1) {
        tr[i].style.display = "";
        } else {
        tr[i].style.display = "none";
        }
    } 
    }
}

function presetsFileSearch() {
    // Declare variables 
    var inputs, filters, table, tr, td, i;
    inputs = [document.getElementById("presetsFileNameInput"), document.getElementById("presetsTFInput"), document.getElementById("presetsTissueInput")];
    filters = [inputs[0].value.toUpperCase(), inputs[1].value.toUpperCase(), inputs[2].value.toUpperCase()];
    table = document.getElementById("presetsFileTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
        tds = tr[i].getElementsByTagName("td");
        if (tds.length > 0) {

            // for (j=0; j < filters.length; j++) {
            // 	console.log(tds[j].innerHTML.toUpperCase().indexOf(filter[j]))
            // }

            if ((tds[0].innerHTML.toUpperCase().indexOf(filters[0]) > -1) 
            && (tds[1].innerHTML.toUpperCase().indexOf(filters[1]) > -1) 
            && (tds[2].innerHTML.toUpperCase().indexOf(filters[2]) > -1)) 
            {
                tr[i].style.display = "";
            } 
            else {
                tr[i].style.display = "none";
            }
        } 
    }
}

function peakFileSearch() {
    // Declare variables 
    var inputs, filters, table, tr, td, i;
    inputs = [document.getElementById("peakFileNameInput"), document.getElementById("peakTFInput"), document.getElementById("peakTissueInput")];
    filters = [inputs[0].value.toUpperCase(), inputs[1].value.toUpperCase(), inputs[2].value.toUpperCase()];
    table = document.getElementById("peakFileTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
        tds = tr[i].getElementsByTagName("td");
        if (tds.length > 0) {

            // for (j=0; j < filters.length; j++) {
            // 	console.log(tds[j].innerHTML.toUpperCase().indexOf(filter[j]))
            // }

            if ((tds[0].innerHTML.toUpperCase().indexOf(filters[0]) > -1) 
            && (tds[1].innerHTML.toUpperCase().indexOf(filters[1]) > -1) 
            && (tds[2].innerHTML.toUpperCase().indexOf(filters[2]) > -1)) 
            {
                tr[i].style.display = "";
            } 
            else {
                tr[i].style.display = "none";
            }
        } 
    }
}

function motifFileSearch() {
    // Declare variables 
    var inputs, filters, table, tr, td, i;
    inputs = [document.getElementById("motifFileNameInput"), document.getElementById("motifTFInput"), document.getElementById("motifTissueInput")];
    filters = [inputs[0].value.toUpperCase(), inputs[1].value.toUpperCase(), inputs[2].value.toUpperCase()];
    table = document.getElementById("motifFileTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
        tds = tr[i].getElementsByTagName("td");
        if (tds.length > 0) {

            // for (j=0; j < filters.length; j++) {
            // 	console.log(tds[j].innerHTML.toUpperCase().indexOf(filter[j]))
            // }

            if ((tds[0].innerHTML.toUpperCase().indexOf(filters[0]) > -1) 
            && (tds[1].innerHTML.toUpperCase().indexOf(filters[1]) > -1) 
            && (tds[2].innerHTML.toUpperCase().indexOf(filters[2]) > -1)) 
            {
                tr[i].style.display = "";
            } 
            else {
                tr[i].style.display = "none";
            }
        } 
    }
}

function motifFileSearch() {
    // Declare variables 
    var inputs, filters, table, tr, td, i;
    inputs = [document.getElementById("enhancerFileNameInput"), document.getElementById("enhancerTissueInput")];
    filters = [inputs[0].value.toUpperCase(), inputs[1].value.toUpperCase()];
    table = document.getElementById("enhancerFileTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
        tds = tr[i].getElementsByTagName("td");
        if (tds.length > 0) {

            // for (j=0; j < filters.length; j++) {
            // 	console.log(tds[j].innerHTML.toUpperCase().indexOf(filter[j]))
            // }

            if ((tds[0].innerHTML.toUpperCase().indexOf(filters[0]) > -1) 
            && (tds[1].innerHTML.toUpperCase().indexOf(filters[1]) > -1)) 
            {
                tr[i].style.display = "";
            } 
            else {
                tr[i].style.display = "none";
            }
        } 
    }
}

function pwmFileSearch() {
    // Declare variables 
    var inputs, filters, table, tr, td, i;
    inputs = [document.getElementById("pwmFileNameInput"), document.getElementById("pwmTFInput")];
    filters = [inputs[0].value.toUpperCase(), inputs[1].value.toUpperCase()];
    table = document.getElementById("pwmFileTable");
    tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (i = 0; i < tr.length; i++) {
        tds = tr[i].getElementsByTagName("td");
        if (tds.length > 0) {

            // for (j=0; j < filters.length; j++) {
            // 	console.log(tds[j].innerHTML.toUpperCase().indexOf(filter[j]))
            // }

            if ((tds[0].innerHTML.toUpperCase().indexOf(filters[0]) > -1) 
            && (tds[1].innerHTML.toUpperCase().indexOf(filters[1]) > -1)) 
            {
                tr[i].style.display = "";
            } 
            else {
                tr[i].style.display = "none";
            }
        } 
    }
}