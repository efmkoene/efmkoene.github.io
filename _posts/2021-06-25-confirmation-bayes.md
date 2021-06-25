---
layout: post
title: "Confirmation Bayes"
subtitle: ""
tags: [psychology,science,short]
---

### Confirmation bias
[Confirmation bias](https://en.wikipedia.org/wiki/Confirmation_bias) is often painted as a bad thing. People supposedly select information (or news sources, ...) that confirm the views they already have, or avoid information that contradicts what they already think. It is considered a bad thing, leading to overconfidence, erroneous decisionmaking in the face of reality, belief in conspiracy theories due to filter bubble effects, and this list goes on. Our ideas would be like babies in our head which we take great care for and protect from conflict.

Its role in science is, in my view, quite interesting. Or, more precisely, *is confirmation bias really a bad thing, or is something else at play?*

### The shame of scientists
An interesting story in this regard is one told by physicist Richard Feynmann:   
> *It's interesting to look at the history of measurements of the charge of an electron, after Millikan. If you plot them as a function of time, you find that one is a little bit bigger than Millikan's, and the next one's a little bit bigger than that, and the next one's a little bit bigger than that, until finally they settle down to a number which is higher.  
Why didn't they discover the new number was higher right away? It's a thing that scientists are ashamed of—this history—because it's apparent that people did things like this: When they got a number that was too high above Millikan's, they thought something must be wrong—and they would look for and find a reason why something might be wrong. When they got a number close to Millikan's value they didn't look so hard. And so they eliminated the numbers that were too far off, and did other things like that...*

Someone actually [tracked down these numbers](https://hsm.stackexchange.com/questions/264/timeline-of-measurements-of-the-electrons-charge) to, indeed, confirm the over-all point raised by Richard Feynmann: you can track the experimentally obtained values for the charge of an electron over time, and they form an interesting line that slowly rises from the wrong value towards the currently accepted value.

![charge of electron](https://i.stack.imgur.com/WtmUj.png)

You can see that Feynmann is describing confirmation bias here: when scientists find something that contradicts their knowledge, they start hunting down problems. Contrarily, when scientists find something that confirms previous beliefs, they don't scrutinize the results so much. But I wouldn't agree with Feynmann that this is something to be ashamed of for various reasons:

1. The entire pursuit of science is driven by the assumption that prior research was done in good faith, and that published results have some *credibility* to them. If scientists would always have to assume that previously found results are incorrect, then how can you develop any theory? Science means standing on the shoulders of giants. Trusting that is not something to be ashamed of.
2. Imagine someone comes around and says their measurement device says that the speed of light is actually twice the currently accepted value. Would it be so wrong for a physicist to simply ignore such evidence or point out that there is probably something wrong with their device, as the finding is just not *credible*? In the words of [Carl Sagan](https://en.wikipedia.org/wiki/Sagan_standard), ''extraordinary claims require extraordinary evidence'', and that is not such a bad standard to have with regards to new theories, is it? A pet-peeve of mine is the [attention Alfred Wegener typically receives for having come up with plate tectonics](https://www.discovermagazine.com/planet-earth/continental-drift-a-revolutionary-theory-that-was-once-considered) while being largely ignored by the established science. The truth is that Alfred Wegener could not provide any appropriate mechanism for plate tectonics. So, in my view, his extraordinary claim was rightly ignored, without the extraordinary evidence! Only when evidence started to appear in the 1950s that supported plate tectonics, did the theory become *credible*.
3. Science is not as objective as sometimes portrayed and popularized. A scientist doesn't start researching things "at random", and then finding things out. No, you set out on a journey, expecting to find a certain thing. Whether you find it or not is tested by experiments etcetera, but the point is that you set out expecting a certain outcome! You hope that your method works better, is more robust, explains more, ..., without having the evidence yet! We have favorite theories, hopes that certain theories are correct and that certain are wrong, even if there is no proof for it! The reality is that science doesn't change in the light of new evidence, but rather that science advances on funeral at a time, paraphrasing [Max Planck](https://www.sciencedaily.com/releases/2019/08/190829150642.htm): 
    
    > "A great scientific truth does not triumph by convincing its opponents and making them see the light, but rather because its opponents eventually die, and a new generation grows up that is familiar with it."  

   One could, of course, still object that this is a bad thing. But I don't think it is. An unhealthy and stubborn expectation that you will find what you expect to find (in the case of the researchers, they expect that their new experiment will reproduce a previously found value) is what keeps you going!
   
### Bias or Bayes?
The above realities of science seem to suggest that confirmation bias is baked into science. But I think a closer inspection may reveal an alternative interpretation. Scientific truth is not something that is settled for once and for all by a single experiment. We think that something is true. Evidence that supports this prior assertion won't be scrutinized much, because it confirms what we know, so it's probably correct.  Conversely, evidence that contradicts what we think we know is heavily derided, heavily criticized, or ignored entirely. 

That may seem ignorant. But I don't think it is. In fact, I think it is a perfectly good strategy following [Bayesian inference](https://en.wikipedia.org/wiki/Bayesian_inference). This is a mathematical theorem describes how **beliefs** should be updated in the face of the likelikhood of new evidence to be true, and old beliefs to be wrong.
