---
layout: page
title: Focusing resource
subtitle: ""
---

<style type="text/css" rel="stylesheet">  
.hidden{ display: none; }
  
.ui-progressbar-value {
  transition: width 0.5s;
  -webkit-transition: width 0.5s;
}

.buttoned {
  background-color: #4CAF50; /* Green */
  border: 1px solid green;
  color: white;
  padding: 7px 16px;
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

  <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
  <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
<link rel="stylesheet" href="//code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">
<script>
$( function() {$( "#progressbar" ).progressbar({value: 10});});
  </script>




I came across *Inner Relationship Focusing* a year ago. I thought it was just weird new wavy bullshit, but I've found it to be actually rather helpful.

My understanding of the method is that you create a space in your body to meet yourself. I know it sounds stupid. But hear me out. You create a space to get a sense on how you're feeling, deep down. The thing is -- your body can't talk. You'll have to interrogate yourself, by asking a lot of "yes" or "no" questions. The body responds with physical sensations every time that you get closer to what it means. This to me shows that there is something to it -- it is beyond just fantasizing, because particular feelings resonate while others don't. The process is not quick (about 15 minutes, maybe), but it is effective. I've used it to figure things out about my responses to events, or figure out how I "really" feel about things.

<h2> Zoom in on how you're feeling </h2>
<!--Timer portion-->
<div id="Q0" class="">
  <h4> 1. Take 30 seconds to just get a sense of yourself. Breath calmly and deeply. How do you feel, in general?</h4>
  <input id="timer1" type="button" value="Start timer..." />
<div id='countdown1'></div>
</div>

<!--Question one-->
<br clear="all" />
<div id="Q1" class="hidden">
  <h4> 2. Do you notice a weight to your body? Or a springiness? </h4>
  <button class="option1 buttoned">Feeling light</button> 
  <button class="option1 buttoned">Intermediate</button>
  <button class="option1 buttoned">Feeling heavy</button>
</div>

<!--Question two-->
<br clear="all" />
<div id="Q2" class="hidden">
  <h4> 3. Where in the body do I sense something taking my attention (e.g., tightness, pressure, a knocking feeling, ...). Where is your body trying to get to speak to you?</h4>
  <button class="option2 buttoned">My chest</button> 
  <button class="option2 buttoned">My belly</button>
  <button class="option2 buttoned">Elsewhere</button>
</div>

<!--Question three; just a timer...-->
<br clear="all" />
<div id="Q3" class="hidden">
  <h4> 4. Try to get a handle on the unclear bodily sensation. Is it a pressure, something moving, something hot/cold, ... . Something vague suffices for now. Try to sense it wholly. Tell yourself "I'm sensing ..., and I am saying hello to that."</h4>
    <div id='countdown2'></div>
</div>

<!--Question four...-->
<br clear="all" />
<div id="Q4" class="hidden">
  <h4> 5. Get a more sense of the emotional quality of the feeling. You do this by testing with yourself. For example, "I'm sensing pain", or "I'm sensing frustration". See what words bubble up. Try to capture the whole feeling in a single qualitative word. Your body will respond if the word resonates with its feeling.</h4>
  <p><img src="https://gritx.org/skills-studio/uploads/exerciseimage/emotion_wheel8.jpg"></p>
  <button class="option4 buttoned">Done.</button>
</div>

<!--Question five...-->
<br clear="all" />
<div id="Q5" class="hidden">
  <h4> 6. Now create a relationship with this felt quality. Say, you'd meet this person in a bar, what kind of bar would it be? Would you be having intense conversations, or would you both sit quietly by? Would you be drinking or eating? Again, use what intuition bubbles up. Check by seeing what resonates. Take a moment again, to imagine this.</h4>
  <button class="option5 buttoned">Done.</button>
</div>

<!--Question six...-->
<br clear="all" />
<div id="Q6" class="hidden">
  <h4> 7. Ask yourself, "What would it feel like if it was all OK?". "What is in the way of that?". See what bubbles up.</h4>
  <div id='countdown3'></div>
</div>

<!--End...-->
<br clear="all" />
<div id="Q7" class="hidden">
  <h4> 8. The end! Thank your body for speaking with you. Remember, you don't have to agree with your body, you just have to acknowledge what it said to you.</h4>
</div>


OK? 

<div id="progressbar"></div>

Ok......


<script> 
$('button').click(function() {
    $(this).css('background-color', 'black');
});
  

var Q0 = document.getElementById('Q0');
var Q1 = document.getElementById('Q1');
var Q2 = document.getElementById('Q2');
var Q3 = document.getElementById('Q3');
var Q4 = document.getElementById('Q4');
var Q5 = document.getElementById('Q5');
var Q6 = document.getElementById('Q6');
var Q7 = document.getElementById('Q7');
var btn1 = document.getElementsByClassName('option1');
var btn2 = document.getElementsByClassName('option2');
var btn4 = document.getElementsByClassName('option4');
var btn5 = document.getElementsByClassName('option5');


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

for(var i=0; i<btn4.length; i++){
    btn4[i].addEventListener("click", function(){ 
  Q5.className = ''; 
})
}
 
for(var i=0; i<btn5.length; i++){
    btn5[i].addEventListener("click", function(){ 
  Q6.className = ''; 
  countdown('countdown3', 120,'Q7');
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
