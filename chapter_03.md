---
# Page settings
layout: default
keywords:
comments: true
image: alaska.jpg
# Hero section
title:  Chapter III
subtitle: >  
    First-principles prediction of the information processing capacity of a simple genetic circuit

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        Living organisms are constantly sensing intra and extracellular cues and
        responding accordingly. The quality and precision of such responses can mean
        the difference between surviving or not certain challenges; therefore, there
        is a constant selection pressure for cells to gather enough information from
        any stimulus in order to build an adequate response. In this context, the
        information that cells can obtain has a precise mathematical definition
        measured-—just as in computers—-in bits.

        In this chapter our goal is to predict how many bits of information can a
        cell harboring a simple genetic circuit process. To do so, I write down a
        theoretical model to predict the full distribution of gene expression based on
        the physics of this molecular process. I calibrate our model with previous
        information in order to perform parameter-free predictions. To test the model, I
        compare the predictions with experimental single-cell gene expression data
        finding great agreement.
# Page navigation
page_nav:
    prev:
        content: Chapter 2
        url: chapter_02
    next:
        content: Chapter 4
        url: chapter_04
prefix: chapter_03
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_minimal_model
    - section_05_distribution_moments
    - section_06_cell_cycle
    - section_07_maxent
    - section_08_channel_capacity
    - section_09_discussion
    - section_10_methods
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
