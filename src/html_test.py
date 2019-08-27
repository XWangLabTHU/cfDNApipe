from markdown.inlinepatterns import Pattern
from markdown.util import etree

class MultiPattern(Pattern):
    def handleMatch(self, m):
        if m.group(2) == '**':
            # Bold
            tag = 'strong'
        elif m.group(2) == '//':
            # Italics
            tag = 'em'
        elif m.group(2) == '__':
            # Underline
            tag = 'ins'
        elif m.group(2) == '--':
            # Strike
            tag = 'del'
        # Create the Element
        el = etree.Element(tag)
        el.text = m.group(3)
        return el