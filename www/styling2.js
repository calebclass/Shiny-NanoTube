var btn1 = document.getElementById("jobSetup");
var btn2 = document.getElementById("QCResults");
var btn3 = document.getElementById("analyzeResults");
var btn4 = document.getElementById("geneSetAnalysis");


btn1.onclick = function(){
  var targetDiv = document.getElementById("btn-group1");
  if (targetDiv.style.display !== "none"){
    targetDiv.style.display = "none";
  } else {
    targetDiv.style.display = "block";
  }
};

btn2.onclick = function(){
  var targetDiv = document.getElementById("btn-group2");
  if (targetDiv.style.display !== "none"){
    targetDiv.style.display = "none";
  } else {
    targetDiv.style.display = "block";
  }
};


btn3.onclick = function(){
  var targetDiv = document.getElementById("btn-group3");
  if (targetDiv.style.display !== "none"){
    targetDiv.style.display = "none";
  } else {
    targetDiv.style.display = "block";
  }
};


btn4.onclick = function(){
  var targetDiv = document.getElementById("btn-group4");
  if (targetDiv.style.display !== "none"){
    targetDiv.style.display = "none";
  } else {
    targetDiv.style.display = "block";
  }
};


