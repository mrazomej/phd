---
# Page settings
layout: default
keywords:
comments: true
image: plants.jpg
# Hero section
title:  Chapter I
subtitle: >  
    Introduction

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        To fully appreciate the beauty of science, many times, it is necessary
        not only to superficially understand concepts behind the natural 
        phenomena being studied, but one needs to dig deeper into the technical
        details that lead us to a deeper understanding of the way the universe 
        works. This is the difference between "knowing about something" and 
        actually "knowing something."

        Not a single piece of work is ever complete, and not a single source of
        information is ever sufficient to understand everything fully. Let alone
        this thesis. Nevertheless, with a few concepts at hand, it is possible 
        to appreciate better the work developed here. In this chapter, I 
        introduce such concepts. I present the conceptual, physical, and 
        mathematical tools I consider as a prerequisite to work through this
        text. From the basics of how to model gene expression in bacteria to
        the ideas from Information Theory, this is an attempt to get the reader
        on the same page, hopefully making the experience of reading the rest
        of the text a much less daunting task. 
         
# Page navigation
page_nav:
    prev:
        content: Abstract
        url: abstract
    next:
        content: Chapter 2
        url: chapter_02
prefix: chapter_01
contents:
    - section_00_title
    - section_01_introduction
    - section_02_gene_reg
    - section_03_information
---

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
