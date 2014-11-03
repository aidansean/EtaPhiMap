from project_module import project_object, image_object, link_object, challenge_object

p = project_object('EtaPhiMap', 'eta-phi map')
p.domain = 'http://www.aidansean.com/'
p.path = 'EtaPhiMap'
p.preview_image_ = image_object('http://placekitten.com.s3.amazonaws.com/homepage-samples/408/287.jpg', 408, 287)
p.github_repo_name = 'EtaPhiMap'
p.mathjax = False
p.introduction = 'This project is intended to make a map from \((x,y)\) space to \((\eta,\phi)\) space to make analysis in the CMS endcaps easier and more intuitive.'
p.overview = '''This script arranges "superclusters" around the endcap and shows what a cone of constant \(\Delta R\) looks like.'''
