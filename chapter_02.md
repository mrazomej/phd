---
# Page settings
layout: default
keywords:
comments: true
image: paihia.jpg
# Hero section
title:  Chapter II
subtitle: >  
    A Predictive Theory of Allosteric Induction

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        Despite lacking a nervous system, single bacterial cells are capable of
        making decisions given signals from their surroundings. How can individual
        molecules sense and transmit these signals? The answer comes from one of the
        crowning scientific achievements of the past century: allostery. Simply
        stated, allostery is the property of certain macromolecules to exist in
        multiple conformations with different properties. For example, transcription
        factors —proteins that control gene expression— can be active (able to bind
        the DNA) or inactive (unable to bind DNA) depending on the concentration of
        a signaling molecule.
        
        In this chapter, I write down a theoretical model that predicts
        the expression level of a gene regulated by an allosteric transcription factor.
        We then test the model experimentally, showing that the model was able to
        predict how changes to the regulation of the gene translate to changes in the
        cellular response.
# Page navigation
page_nav:
    prev:
        content: Chapter 1
        url: chapter_01
    next:
        content: Chapter 3
        url: chapter_03
prefix: chapter_02
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_model
    - section_05_results
    - section_06_discussion
    - section_07_methods
---

**Published as ...**
<hr/>
{% if page.contents %}
{% for val in page.contents %}
{% if jekyll.environment == production %}
{% include_relative {{site.doks.baseurl}}src/{{page.prefix}}/{{val}}.md %}
{% else %}
{% include_relative src/{{page.prefix}}/{{val}}.md %}
{% endif %}
{% endfor %}
{% endif %}

## References
