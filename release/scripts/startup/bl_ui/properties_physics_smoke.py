# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>
import bpy
import os
from copy import deepcopy
from bpy.types import Panel

from bl_ui.properties_physics_common import (
        point_cache_ui,
        effector_weights_ui,
        )


class PhysicButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "physics"

    @classmethod
    def poll(cls, context):
        ob = context.object
        rd = context.scene.render
        return (ob and ob.type == 'MESH') and (not rd.use_game_engine) and (context.smoke)


class PHYSICS_PT_smoke(PhysicButtonsPanel, Panel):
    bl_label = "Smoke"

    def draw(self, context):
        layout = self.layout

        md = context.smoke
        ob = context.object

        layout.prop(md, "smoke_type", expand=True)

        if md.smoke_type == 'DOMAIN':
            domain = md.domain_settings

            split = layout.split()

            split.enabled = not domain.point_cache.is_baked

            col = split.column()
            col.label(text="Resolution:")
            col.prop(domain, "resolution_max", text="Divisions")
            col.label(text="Time:")
            col.prop(domain, "time_scale", text="Scale")
            col.label(text="Border Collisions:")
            col.prop(domain, "collision_extents", text="")

            col = split.column()
            col.label(text="Behavior:")
            col.prop(domain, "alpha")
            col.prop(domain, "beta", text="Temp. Diff.")
            col.prop(domain, "vorticity")
            col.prop(domain, "use_dissolve_smoke", text="Dissolve")
            sub = col.column()
            sub.active = domain.use_dissolve_smoke
            sub.prop(domain, "dissolve_speed", text="Time")
            sub.prop(domain, "use_dissolve_smoke_log", text="Slow")

        elif md.smoke_type == 'FLOW':

            flow = md.flow_settings

            layout.prop(flow, "smoke_flow_type", expand=False)

            if flow.smoke_flow_type != 'OUTFLOW':
                split = layout.split()
                col = split.column()
                col.label(text="Flow Source:")
                col.prop(flow, "smoke_flow_source", expand=False, text="")
                if flow.smoke_flow_source == 'PARTICLES':
                    col.label(text="Particle System:")
                    col.prop_search(flow, "particle_system", ob, "particle_systems", text="")
                    col.prop(flow, "use_particle_size", text="Set Size")
                    sub = col.column()
                    sub.active = flow.use_particle_size
                    sub.prop(flow, "particle_size")
                else:
                    col.prop(flow, "surface_distance")
                    col.prop(flow, "volume_density")

                sub = col.column(align=True)
                sub.prop(flow, "use_initial_velocity")

                sub = sub.column()
                sub.active = flow.use_initial_velocity
                sub.prop(flow, "velocity_factor")
                if flow.smoke_flow_source == 'MESH':
                    sub.prop(flow, "velocity_normal")
                    #sub.prop(flow, "velocity_random")

                sub = split.column()
                sub.label(text="Initial Values:")
                sub.prop(flow, "use_absolute")
                if flow.smoke_flow_type in {'SMOKE', 'BOTH'}:
                    sub.prop(flow, "density")
                    sub.prop(flow, "temperature")
                    sub.prop(flow, "smoke_color")
                if flow.smoke_flow_type in {'FIRE', 'BOTH'}:
                    sub.prop(flow, "fuel_amount")
                sub.label(text="Sampling:")
                sub.prop(flow, "subframes")

        elif md.smoke_type == 'COLLISION':
            coll = md.coll_settings

            split = layout.split()

            col = split.column()
            col.prop(coll, "collision_type")


class PHYSICS_PT_smoke_flow_advanced(PhysicButtonsPanel, Panel):
    bl_label = "Smoke Flow Advanced"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        return md and (md.smoke_type == 'FLOW') and (md.flow_settings.smoke_flow_source == 'MESH')

    def draw(self, context):
        layout = self.layout
        ob = context.object
        flow = context.smoke.flow_settings

        split = layout.split()
        col = split.column()

        col.prop(flow, "use_texture")
        sub = col.column()
        sub.active = flow.use_texture
        sub.prop(flow, "noise_texture", text="")
        sub.label(text="Mapping:")
        sub.prop(flow, "texture_map_type", expand=False, text="")
        if flow.texture_map_type == 'UV':
            sub.prop_search(flow, "uv_layer", ob.data, "uv_textures", text="")
        if flow.texture_map_type == 'AUTO':
            sub.prop(flow, "texture_size")
        sub.prop(flow, "texture_offset")

        col = split.column()
        col.label(text="Vertex Group:")
        col.prop_search(flow, "density_vertex_group", ob, "vertex_groups", text="")


class PHYSICS_PT_smoke_fire(PhysicButtonsPanel, Panel):
    bl_label = "Smoke Flames"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        return md and (md.smoke_type == 'DOMAIN')

    def draw(self, context):
        layout = self.layout
        domain = context.smoke.domain_settings

        split = layout.split()
        split.enabled = not domain.point_cache.is_baked

        col = split.column(align=True)
        col.label(text="Reaction:")
        col.prop(domain, "burning_rate")
        col.prop(domain, "flame_smoke")
        col.prop(domain, "flame_vorticity")

        col = split.column(align=True)
        col.label(text="Temperatures:")
        col.prop(domain, "flame_ignition")
        col.prop(domain, "flame_max_temp")
        col.prop(domain, "flame_smoke_color")


class PHYSICS_PT_smoke_adaptive_domain(PhysicButtonsPanel, Panel):
    bl_label = "Smoke Adaptive Domain"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        return md and (md.smoke_type == 'DOMAIN')

    def draw_header(self, context):
        md = context.smoke.domain_settings

        self.layout.prop(md, "use_adaptive_domain", text="")
		
    def draw(self, context):
        layout = self.layout

        domain = context.smoke.domain_settings
        layout.active = domain.use_adaptive_domain
		
        split = layout.split()
        split.enabled = (not domain.point_cache.is_baked)

        col = split.column(align=True)
        col.label(text="Resolution:")
        col.prop(domain, "additional_res")
        col.prop(domain, "adapt_margin")

        col = split.column(align=True)
        col.label(text="Advanced:")
        col.prop(domain, "adapt_threshold")


class PHYSICS_PT_smoke_highres(PhysicButtonsPanel, Panel):
    bl_label = "Smoke High Resolution"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        rd = context.scene.render
        return md and (md.smoke_type == 'DOMAIN') and (not rd.use_game_engine)

    def draw_header(self, context):
        md = context.smoke.domain_settings

        self.layout.prop(md, "use_high_resolution", text="")

    def draw(self, context):
        layout = self.layout

        md = context.smoke.domain_settings

        layout.active = md.use_high_resolution

        split = layout.split()
        split.enabled = not md.point_cache.is_baked

        col = split.column()
        col.label(text="Resolution:")
        col.prop(md, "amplify", text="Divisions")
        col.label(text="Flow Sampling:")
        col.row().prop(md, "highres_sampling", text="")

        col = split.column()
        col.label(text="Noise Method:")
        col.row().prop(md, "noise_type", text="")
        col.prop(md, "strength")

        layout.prop(md, "show_high_resolution")


class PHYSICS_PT_smoke_groups(PhysicButtonsPanel, Panel):
    bl_label = "Smoke Groups"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        rd = context.scene.render
        return md and (md.smoke_type == 'DOMAIN') and (not rd.use_game_engine)

    def draw(self, context):
        layout = self.layout
        domain = context.smoke.domain_settings

        split = layout.split()

        col = split.column()
        col.label(text="Flow Group:")
        col.prop(domain, "fluid_group", text="")

        #col.label(text="Effector Group:")
        #col.prop(domain, "effector_group", text="")

        col = split.column()
        col.label(text="Collision Group:")
        col.prop(domain, "collision_group", text="")


class PHYSICS_PT_smoke_cache(PhysicButtonsPanel, Panel):
    bl_label = "Smoke Cache"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        rd = context.scene.render
        return md and (md.smoke_type == 'DOMAIN') and (not rd.use_game_engine)

    def draw(self, context):
        layout = self.layout

        md = context.smoke.domain_settings
        cache = md.point_cache

        layout.label(text="Compression:")
        layout.prop(md, "point_cache_compress_type", expand=True)

        point_cache_ui(self, context, cache, (cache.is_baked is False), 'SMOKE')


class PHYSICS_PT_smoke_field_weights(PhysicButtonsPanel, Panel):
    bl_label = "Smoke Field Weights"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.smoke
        rd = context.scene.render
        return md and (md.smoke_type == 'DOMAIN') and (not rd.use_game_engine)

    def draw(self, context):
        domain = context.smoke.domain_settings
        effector_weights_ui(self, context, domain.effector_weights, 'SMOKE')

class OBJECT_OT_RunMantaButton(bpy.types.Operator):
    bl_idname = "manta_export_scene.button"
    bl_label = "Create Python Script and mesh files"
    
    def execute(self, context):
        def silent_remove(filename):
            if os.path.exists(filename):
                os.remove(filename)

        #need these methods to account for rotated objects
        def transform_objgroup(obj_list, domain_obj):
            old_scale = deepcopy(domain_obj.scale)
            old_loc = deepcopy(domain_obj.location)
            #link all objects to new reference- domain
            domain_obj.scale = (1,1,1)
            domain_obj.location = (0,0,0)
            for obj in obj_list:
                obj.select = True
                obj.location[0] -= old_loc[0]
                obj.location[1] -= old_loc[1]
                obj.location[2] -= old_loc[2]
                obj.constraints.new('CHILD_OF')
                obj.constraints.active.target = domain_obj
            #scale domain down
            domain_obj.scale[0] /= old_scale[0]
            domain_obj.scale[1] /= old_scale[1]
            domain_obj.scale[2] /= old_scale[2]
            return old_scale, old_loc
            
        def transform_objgroup_back(obj_list, domain_obj, old_data):
            old_scale, old_loc = old_data
            domain_obj.scale[0] =  old_scale[0]
            domain_obj.scale[1] =  old_scale[1]
            domain_obj.scale[2] =  old_scale[2]
            domain_obj.location[0] =  old_loc[0]
            domain_obj.location[1] =  old_loc[1]
            domain_obj.location[2] =  old_loc[2]
            #remove used constraint and deselect objects
            for obj in obj_list:
                obj.select = False
                obj.constraints.remove(obj.constraints.active)
                obj.location[0] += old_loc[0]
                obj.location[1] += old_loc[1]
                obj.location[2] += old_loc[2]
        
        coll_objs = []
        flow_objs = []
        selected_before = []
        domain = None
        #getting smoke objects
        for scene in bpy.data.scenes:
            for ob in scene.objects:
                for modifier in ob.modifiers:
                    if modifier.type == 'SMOKE':
                        if modifier.smoke_type == 'COLLISION':
                            coll_objs.append(ob)
                        elif modifier.smoke_type == 'FLOW':
                            flow_objs.append(ob)
                        elif modifier.smoke_type == 'DOMAIN' and ob.select:
                            domain = ob
                if ob.select:
                    selected_before.append(ob)
                    ob.select = False
        manta_filepath = os.path.dirname(context.smoke.domain_settings.manta_filepath)
        silent_remove(os.path.join(manta_filepath, "manta_coll.obj"))
        silent_remove(os.path.join(manta_filepath, "manta_flow.obj"))
        #exporting here
        if coll_objs: 
            old_data = transform_objgroup(coll_objs, domain)
            bpy.ops.export_scene.obj(filepath = os.path.join(manta_filepath, "manta_coll.obj"), axis_forward='Y', axis_up='Z', use_selection = True, use_normals = True, use_materials = False, use_triangles = True, group_by_object = True, use_nurbs=True, check_existing= False)
            transform_objgroup_back(coll_objs,domain,old_data)
        if flow_objs:
            old_data = transform_objgroup(flow_objs, domain)
            bpy.ops.export_scene.obj(filepath = os.path.join(manta_filepath, "manta_flow.obj"), axis_forward='Y', axis_up='Z', use_selection = True, use_normals = True, use_materials = False, use_triangles = True, group_by_object = True, use_nurbs=True, check_existing= False)
            transform_objgroup_back(flow_objs,domain,old_data)
        for ob in selected_before:
            ob.select = True
        # ds = domain.modifiers['Smoke'].domain_settings
        # if (!global manta_solver_res_switched) and ds.manta_solver_res == 2:
        #     #resize domain s.th. Y-axis dim corresponds to 1
        #     scale_fac = ds.resolution_max / max(domain.scale[0],domain.scale[1],domain.scale[2])
        #     domain.scale[1] /= scale_fac
        #     global manta_solver_res_switched = True
        bpy.ops.manta.make_file()
        bpy.ops.manta.sim_step()
        return{'FINISHED'}

class OBJECT_OT_StopMantaButton(bpy.types.Operator):
    bl_idname = "manta_stop_sim.button"
    bl_label = "Stop Mantaflow Simulation"
    def execute(self, context):
        domain = context.smoke.domain_settings
        #setting manta_sim_frame to "stop" value 
        domain.manta_sim_frame = -1
        return{'FINISHED'}


class PHYSICS_PT_smoke_manta_settings(PhysicButtonsPanel, Panel):
    bl_label = "MantaFlow Settings"
    bl_options = {'DEFAULT_CLOSED'}
    name = bpy.props.StringProperty(name="Test Prop", default="Unknown")
    StringProp = bpy.props.StringProperty(name="manta_status", description="Status Of Simulation", default="Doing Nothing" )
#    filepath = StringProperty(subtype='FILE_PATH',)
    @classmethod
    def poll(cls, context):
        md = context.smoke
        return md and (md.smoke_type == 'DOMAIN')
    
    def draw_header(self, context):
        md = context.smoke.domain_settings
        
    def draw(self, context):
        layout = self.layout
        
        domain = context.smoke.domain_settings
        split = layout.split()
        split.prop(domain, "use_manta_liquid", text="Liquid")
        split.operator("manta_export_scene.button", text="Create Manta Setup")
        split = layout.split()
        split.prop(domain, "manta_filepath")
        split = layout.split()
        col = split.column()
        col.prop(domain, "manta_solver_res", text="Solver Resolution")
        col.prop(domain, "manta_uvs", text="UVs count")
        split = layout.split()
        col = split.column()
        col.label("Noise Settings")
        col.prop(domain, "noise_clamp_neg", text="Clamp Neg")
        col.prop(domain, "noise_clamp_pos", text="Clamp Pos")
        col.prop(domain, "noise_time_anim", text="Time Anim")
        col = split.column()
        col.label("")
        col.prop(domain, "noise_val_scale", text="Scale")
        col.prop(domain, "noise_val_offset", text="Offset")

if __name__ == "__main__":  # only for live edit.
    bpy.utils.register_module(__name__)
