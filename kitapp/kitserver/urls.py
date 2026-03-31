from django.urls import path
from . import views

urlpatterns = [
    path("index",views.index,name = "index"),
    path("register_vistor", views.register_visitor, name="register_vistor"),
    path("register_visitor", views.register_visitor, name="register_visitor"),
    path("track_visit", views.track_visit, name="track_visit"),
    path("initdata",views.InitData, name="initdata"),
    path("assembly",views.Assembly,name = "assembly"),
    path("task_status/<str:taskID>",views.task_status,name="task_status"),
    path("getAssembly/<str:taskID>/<str:name>",views.getAssembly,name="getAssembly"),
    path("gg_assemble",views.gg_assemble,name="gg_assemble"),
    path("getTutorial",views.getTutorial,name="getTutorial"),
    path("getZip",views.getZip,name="getZip"),
]
