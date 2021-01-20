---
layout: page
title: Focusing resource
subtitle: ""
---

<style type="text/css" rel="stylesheet">  
.hidden{ display: none; }

.button {
  background-color: #4CAF50; /* Green */
  border: 1px solid green;
  color: white;
  padding: 15px 32px;
  text-align: center;
  text-decoration: none;
  display: inline-block;
  font-size: 16px;
  cursor: pointer;
  float: left;
}

.base-timer {
  position: relative;
  width: 300px;
  height: 300px;
}

.base-timer__svg {
  transform: scaleX(-1);
}

.base-timer__circle {
  fill: none;
  stroke: none;
}

.base-timer__path-elapsed {
  stroke-width: 7px;
  stroke: grey;
}

.base-timer__path-remaining {
  stroke-width: 7px;
  stroke-linecap: round;
  transform: rotate(90deg);
  transform-origin: center;
  transition: 1s linear all;
  fill-rule: nonzero;
  stroke: currentColor;
}

.base-timer__path-remaining.green {
  color: rgb(65, 184, 131);
}

.base-timer__path-remaining.orange {
  color: orange;
}

.base-timer__path-remaining.red {
  color: red;
}

.base-timer__label {
  position: absolute;
  width: 300px;
  height: 300px;
  top: 0;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 48px;
}
</style>




I came across *Inner Relationship Focusing* a year ago. I thought it was just weird new wavy bullshit, but I've found it to be actually rather helpful.

My understanding of the method is that you create a space in your body to meet yourself. I know it sounds stupid. But hear me out. You create a space to get a sense on how you're feeling, deep down. 

V29.

<h2> Zoom in on how you're feeling </h2>
<!--Timer portion-->
<div id="Q0" class="">
  <h4> Take 30 seconds to just get a sense of yourself. Breath calmly and deeply. How do you feel, in general?</h4>
  <input id="timer1" type="button" value="Start timer..." />
<div id='countdown1'></div>
</div>

<!--Question one-->
<div id="Q1" class="hidden">
  <h4> Do you notice a weight to your body? Or a springiness? </h4>
  <button class="option1 button">Feeling light</button> 
  <button class="option1 button">Intermediate</button>
  <button class="option1 button">Feeling heavy</button>
</div>

<!--Question two-->
<div id="Q2" class="hidden">
  <h4> Where in the body do I sense something taking my attention (e.g., tightness, pressure, a knocking feeling, ...)</h4>
  <button class="option2 button">My chest</button> 
  <button class="option2 button">My belly</button>
  <button class="option2 button">Elsewhere</button>
</div>

<!--Question three, with new timer...-->
<div id="Q3" class="hidden">
  <h4> Welcome the feeling, say "I am sensing ... in me, and I'm saying hello to that feeling"</h4>
  <div id='countdown2'></div>
</div>

<div id="Q4" class="hidden">
  <h4> "Try and characterize the feeling. For example, " </h4>
</div>


<script src="https://code.jquery.com/jquery-3.5.1.slim.min.js" integrity="sha256-4+XzXVhsDmqanXGHaHvgh1gMQKX40OUvDEBTu8JcmNs=" crossorigin="anonymous"></script>
<script>
$('button').click(function() {
    $(this).css('background-color', black);
    $(this).effect( "highlight", {color: #4CAF50}, 3000 );
});

var Q0 = document.getElementById('Q0');
var Q1 = document.getElementById('Q1');
var Q2 = document.getElementById('Q2');
var btn1 = document.getElementsByClassName('option1');
var btn2 = document.getElementsByClassName('option2');

for(var i=0; i<btn1.length; i++){
    btn1[i].addEventListener("click", function(){ 
  Q2.className = ''; 
})
}

for(var i=0; i<btn2.length; i++){
    btn2[i].addEventListener("click", function(){ 
  Q3.className = ''; 
  countdown('countdown2', 20,'Q4');
})
}
  
  
  
  
  
  
  
  
  
  
// Countdown timer stuff  
function countdown(element, seconds, next_class) {
    // Fetch the display element
    seconds = seconds*100;
    var total_time=seconds;
    var el = document.getElementById(element);
    
    var nex = document.getElementById(next_class)
    
    
    

    // Set the timer
    var interval = setInterval(function() {
        if(seconds == 0) {
            clearInterval(interval);
            return;
        }
        if(seconds < (5)*100) {
          nex.className = '';
        }
        el.innerHTML = `
<div class="base-timer";>
  <svg class="base-timer__svg" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">
    <g class="base-timer__circle">
      <circle class="base-timer__path-elapsed" cx="50" cy="50" r="25"></circle>
      <path id="base-timer-path-remaining" stroke-dasharray="`+ (seconds)/total_time*157 + ` 157" class="base-timer__path-remaining green" d="
          M 30, 30
          m -5, 20
          a 25,25 0 1,0 50,0
          a 25,25 0 1,0 -50,0
        "></path>
    </g>
  </svg>
  <span id="base-timer-label" class="base-timer__label">`+ Math.floor(seconds/100+0.5) +`</span>
</div>`;
        
    seconds--;
    }, 10); // Update every 10 ms
}

function fixedcount(element, seconds) {
    // Fetch the display element
    seconds = seconds*100;
    var total_time=seconds;
    var el = document.getElementById(element);
    // Set the timer
        el.innerHTML = `
<div class="base-timer";>
  <svg class="base-timer__svg" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">
    <g class="base-timer__circle">
      <circle class="base-timer__path-elapsed" cx="50" cy="50" r="25"></circle>
      <path id="base-timer-path-remaining" stroke-dasharray="`+ (seconds)/total_time*157 + ` 157" class="base-timer__path-remaining green" d="
          M 30, 30
          m -5, 20
          a 25,25 0 1,0 50,0
          a 25,25 0 1,0 -50,0
        "></path>
    </g>
  </svg>
  <span id="base-timer-label" class="base-timer__label">`+ Math.floor(seconds/100+0.5) +`</span>
</div>`;
}



// Make buttons load timers
var start1 = document.getElementById('timer1');

fixedcount('countdown1', 30) 

start1.onclick = function() {
    countdown('countdown1', 30,'Q1');
}
</script>
