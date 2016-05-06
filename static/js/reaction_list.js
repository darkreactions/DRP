            
function createRxnDiv(item, index) {
    var firstReaction = document.getElementsByClassName("reaction")[0];
    var newReaction = document.createElement('div'); 
    newReaction.setAttribute('class', 'reaction');
    newReaction.innerHTML = 'This will be reaction' + item.reference;
    firstReaction.parentNode.appendChild(newReaction);
}

function fetchRxn(ev) {
    if(!fetchRxn.currentPage) { 
        fetchRxn.currentPage = 2;
    }
    if((window.innerHeight + window.scrollY) > document.body.offsetHeight*0.5 && !fetchRxn.waiting) {
        fetchRxn.waiting=true;
        var request = new XMLHttpRequest();
        request.onreadystatechange = function() {
            if (request.readyState == 4 && request.status == 200) {
                reactions = JSON.parse(request.responseText);
                fetchRxn.currentPage++;
                reactions.reactions.forEach(createRxnDiv);
                fetchRxn.waiting = false;
            }
        }
        request.open("GET", "http://darkreactions.lavoisier.haverford.edu/database.json?page=" + fetchRxn.currentPage.toString());
        request.send();
    }
}
window.onscroll = fetchRxn;
