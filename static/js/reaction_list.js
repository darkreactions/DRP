function fetchRxn(ev) {
    if(!fetchRxn.currentPage) { 
        fetchRxn.currentPage = 2;
    }
    if(((window.innerHeight + window.scrollY) > document.body.offsetHeight*0.5 || (window.innerHeight + window.pageYOffset) > document.body.offsetHeight*0.5) && !fetchRxn.waiting) {
        fetchRxn.waiting=true;
        var request = new XMLHttpRequest();
        request.onreadystatechange = function() {
            if (request.readyState == 4 && request.status == 200) {
                fetchRxn.currentPage++;
                var firstReaction = document.getElementsByClassName("reaction")[0];
                firstReaction.parentElement.innerHTML = firstReaction.parentElement.innerHTML + request.responseText; 
                fetchRxn.waiting = false;
            }
        }
        request.open("GET", "/database.html?reactions_only=1&page=" + fetchRxn.currentPage.toString());
        request.send();
    }
}
window.onscroll = fetchRxn;
