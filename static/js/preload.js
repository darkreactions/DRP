//Refresh screen.
function refreshScreen(rememberPage) {
	alert("1");
	rememberPage = rememberPage !== undefined ? rememberPage : true
	alert("2");
	if (rememberPage) { 
		setPageCookie(); 
	}
	alert("4");
	window.setTimeout('location.reload()', 300);
}
